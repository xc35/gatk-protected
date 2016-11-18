package org.broadinstitute.hellbender.tools.walkers.contest;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.variant.variantcontext.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Estimate cross-sample contamination
 *
 * This tool determine the percent contamination of an input bam by sample, by lane, or in aggregate across all the input reads.
 *
 * <h3>Usage examples</h3>
 * <p>These are example commands that show how to run ContEst for typical use cases. Square brackets ("[ ]")
 * indicate optional arguments. Note that parameter values and/or resources shown here may not be the latest recommended; see the Best Practices documentation for detailed recommendations. </p>
 *
 * <br />
 * <h4>Contamination estimation using the normal BAM for genotyping on-the-fly</h4>
 * <pre>
 *   ./gatk-launch
 *     ContEst \
 *     -R reference.fasta \
 *     -I estimate_contamination.bam \
 *     -exac exac.vcf \
 *     -L populationSites.interval_list \\TODO: do we need this arg and the next?
 *     [-L targets.interval_list] \
 *     -isr INTERSECTION \  //TODO: what's this?
 *     -o output.txt
 * </pre>
 *
 *<h3>Output</h3>
 * A text file containing estimated percent contamination, as well as error bars on this estimate.
 *
 * <h3>Notes</h3>
 * Multiple modes are supported simultaneously, e.g. contamination by sample and readgroup can be computed in the same run.
 */
@CommandLineProgramProperties(
        summary = " Estimate cross-sample contamination\n" +
                " * This tool determine the percent contamination of an input bam by sample, by lane, or in aggregate across all the input reads.",
        oneLineSummary = "Estimate cross-sample contamination",
        programGroup = QCProgramGroup.class)
public final class ContEst extends LocusWalker {
    protected static final Logger logger = LogManager.getLogger(ContEst.class);

    public enum AggregationLevel { SAMPLE, READGROUP, BAM }

    // the population information; the allele frequencies for each position in known populations
    @Argument(shortName = "exac", doc="subsetted exac vcf containing population allele frequencies", optional = false)
    public FeatureInput<VariantContext> exac;

    @Argument(fullName = "evalSampleName", optional = false, doc = "eval sample name as in BAM header")
    public String evalSample;

    @Argument(fullName = "min_qscore", optional = true, doc = "threshold for minimum base quality score")
    public int MIN_QSCORE = 20;

    @Argument(fullName = "min_mapq", optional = true, doc = "threshold for minimum mapping quality score")
    public int MIN_MAPQ = 20;

    @Argument(shortName = "agg", fullName = "aggregations", doc = "set to BAM (default), SAMPLE or READGROUP to produce per-bam, per-sample or per-lane estimates", optional = true)
    private Set<AggregationLevel> aggregations = null;

    @Argument(shortName = "pc", fullName = "precision", doc = "the degree of precision to which the contamination tool should estimate (e.g. the bin size)", optional = true)
    private double precision = 0.01;

    private static final Allele[] ALLELES = {Allele.create((byte) 'A'), Allele.create((byte) 'C'), Allele.create((byte) 'G'), Allele.create((byte) 'T')};

    Collection<String> allSamples;
    Collection<String> allReadGroups;
    boolean aggregateAtBamLevel;
    boolean aggregateAtSampleLevel;
    boolean aggregateAtReadGroupLevel;

    int countGenotypeHomVar = 0;

    // a map of our contamination targets and their cumulativeEstimates
    // key: aggregation entity ("BAM", sample name, or read group ID)
    // value: cumulative ContaminationEstimate
    private Map<String, ContaminationEstimate> cumulativeEstimateByAggregationKey = new HashMap<>();

    @Override
    public void onTraversalStart() {
        final SampleList samplesList = new IndexedSampleList(new ArrayList<>(ReadUtils.getSamplesFromHeader(getHeaderForReads())));
        if (!samplesList.asListOfSamples().contains(evalSample)) {
            throw new UserException.BadInput("BAM header sample names " + samplesList.asListOfSamples() + "does not contain given tumor" + " sample name " + evalSample);
        }

        aggregateAtBamLevel = aggregations == null || aggregations.contains(AggregationLevel.BAM);
        aggregateAtSampleLevel = aggregations.contains(AggregationLevel.SAMPLE);
        aggregateAtReadGroupLevel = aggregations.contains(AggregationLevel.READGROUP);

        allReadGroups = getHeaderForReads().getReadGroups().stream().map(SAMReadGroupRecord::getId).collect(Collectors.toList());
        allSamples = getHeaderForReads().getReadGroups().stream().map(SAMReadGroupRecord::getSample).distinct().collect(Collectors.toList());
    }

    /**
     * our map function, which emits a contamination stats for each of the subgroups (lanes, samples, etc) that we encounter
     *
     * @param ref     the reference information at this position
     * @param context the read context, where we get the alignment data
     * @param features the reference meta data tracker, from which we get the array truth data
     * @return a mapping of our subgroup name to contamination estimate
     */
    @Override
    public void apply(final AlignmentContext context, final ReferenceContext ref, final FeatureContext features) {
        final List<VariantContext> vcs = features.getValues(exac);
        if (vcs.isEmpty()) {
            return;
        }

        final VariantContext popVC = vcs.get(0);
        final Genotype genotype = getGenotype(context, ref, getHeaderForReads());

        // only use homozygous sites
        if (genotype == null || !genotype.isHomVar()) {
            return;
        } else {
            countGenotypeHomVar++;
        }

        // only use non-reference sites
        final byte myBase = genotype.getAllele(0).getBases()[0];


        final ReadPileup pileup = context.getBasePileup().getPileupForSample(evalSample, getHeaderForReads())
                .makeFilteredPileup(pe -> pe.getQual() >= MIN_QSCORE && pe.getMappingQual() >= MIN_MAPQ);


        final Map<String, ContaminationEstimate> resultsForThisSite = new LinkedHashMap<>();
        if (aggregateAtBamLevel) {
            final ContaminationEstimate results = calcStats(pileup, myBase, popVC);
            addSiteEstimate("BAM", results);
        }

        if (aggregateAtReadGroupLevel) {
            for (final String rg : allReadGroups) {
                final ReadPileup rgPileup = pileup.makeFilteredPileup(pe -> pe.getRead().getReadGroup().equals(rg));
                final ContaminationEstimate results = calcStats(rgPileup, myBase, popVC);
                addSiteEstimate(rg, results);
            }
        }

        if (aggregateAtSampleLevel) {
            for (final String sample : allReadGroups) {
                final ReadPileup samplePileup = pileup.getPileupForSample(sample, getHeaderForReads())
                final ContaminationEstimate results = calcStats(samplePileup, myBase, popVC);
                addSiteEstimate(sample, results);
            }
        }

        for (final String aggregationKey : resultsForThisSite.keySet()) {
            final ContaminationEstimate newStats = resultsForThisSite.get(aggregationKey);

            // merge the sets
            if (cumulativeEstimateByAggregationKey.containsKey(aggregationKey)) {
                cumulativeEstimateByAggregationKey.get(aggregationKey).addSiteData(newStats);
            } else {
                cumulativeEstimateByAggregationKey.put(aggregationKey, newStats);
            }
        }

    }

    // siteEstimate can be null, in which case nothing happens
    private void addSiteEstimate(final String aggregationKey, final ContaminationEstimate siteEstimate) {
        if (siteEstimate == null) {
            return;
        }
        if (cumulativeEstimateByAggregationKey.containsKey(aggregationKey)) {
            cumulativeEstimateByAggregationKey.get(aggregationKey).addSiteData(siteEstimate);
        } else {
            cumulativeEstimateByAggregationKey.put(aggregationKey, siteEstimate);
        }
    }

    private Genotype getGenotype(final AlignmentContext context, final ReferenceContext referenceContext, final SAMFileHeader header) {
        //TODO: or the whole pileup if evalSample is not given?
        final ReadPileup pileup = context.getBasePileup().getPileupForSample(evalSample, header);
        if (pileup == null || pileup.isEmpty()) {
            return null;
        }

        //TODO: this is a really dumb way of calling genotypes.  Replace with a posterior probability.  Then it won't be
        //TODO: necessary to use a minimum depth, which I have already gotten rid of in anticipation
        final double MIN_GENOTYPE_RATIO = 0.8;
        final int[] bases = pileup.getBaseCounts();
        final int maxAlleleIndex = GATKProtectedMathUtils.maxIndex(bases);
        final int numReads = (int) MathUtils.sum(bases);
        final String refBase = String.valueOf((char)referenceContext.getBase());
        if (bases[maxAlleleIndex] / (double) numReads >= MIN_GENOTYPE_RATIO && !refBase.equals(ALLELES[maxAlleleIndex].getBaseString())) {
            return new GenotypeBuilder(evalSample, Collections.singletonList(ALLELES[maxAlleleIndex])).make();
        }
        return null;
    }

    private static final class PopulationFrequencyInfo {
        private byte majorAllele;
        private byte minorAllele;
        private double minorAlleleFrequency;

        private PopulationFrequencyInfo(final byte majorAllele, final byte minorAllele, final double minorAlleleFrequency) {
            this.majorAllele = majorAllele;
            this.minorAllele = minorAllele;
            this.minorAlleleFrequency = minorAlleleFrequency;
        }

        public byte getMajorAllele() {
            return majorAllele;
        }

        public byte getMinorAllele() {
            return minorAllele;
        }

        public double getMinorAlleleFrequency() {
            return minorAlleleFrequency;
        }
    }

    private static PopulationFrequencyInfo parseExACInfo(final VariantContext variantContext) {
        // get the ref and alt allele from an ExAC (subsetted) vcf VariantContext
        // for simplicity, we should subset to alleles that are uncommon (maf < some threshold like 10%) in all populations
        // that way we can simply take the overall human maf

        //TODO: fill this in!!!
        return new PopulationFrequencyInfo(majorAllele, minorAllele, maf);
    }


    /**
     * Calculate the contamination values per division, be it lane, meta, sample, etc
     * @param pileup the pileup
     * @param myAllele the allele we have (our hom var genotype allele)
     * @param exacVC the population variant context from hapmap
     * @return a mapping of each target population to their estimated contamination
     */
    private ContaminationEstimate calcStats(final ReadPileup pileup, final byte myAllele, final VariantContext exacVC) {
            final PopulationFrequencyInfo info = parseExACInfo(exacVC);
            final double alleleFreq = info.getMinorAlleleFrequency();

            final long majorCounts = counter.get(Nucleotide.valueOf(info.getMajorAllele()));
            final long minorCounts = counter.get(Nucleotide.valueOf(info.getMinorAllele()));

            // only use sites where this is the minor allele
            if (myAllele == info.minorAllele) {
                return new ContaminationEstimate(minorCounts, majorCounts, counter, precision, alleleFreq, pileup, info.getMinorAllele(), info.getMajorAllele());
            } else {
                return null;
            }
    }

    @Override
    public Object onTraversalSuccess() {
        outputReport();
        logger.info("Homozygous variant sites: " + countGenotypeHomVar);
        return "SUCCESS";
    }

    public void outputReport() {

        //TODO: logger.info isn't correct -- use a TableUtils method
        logger.info("name\tpopulation\tpopulation_fit\tcontamination\tconfidence_interval_95_width\tconfidence_interval_95_low\tconfidence_interval_95_high\tsites");

        //TODO: very confused why cumulativeEstimates = entry.getValue isn't used
        // TODO: and why have both storedData and cumulativeEstimates???
        for (final Map.Entry<String, ContaminationEstimate> entry : cumulativeEstimateByAggregationKey.entrySet()) {
            final ContaminationEstimate data = entry.getValue();
            final String aggregationKey = entry.getKey();

            final String pm = "%3." + Math.round(Math.log10(1/precision)) +"f";

            final int bins = data.getBins().length;

            //TODO: sort before output?

            final double[][] matrix = new double[bins][data.size()];

            for (int i = 0; i<bins; i++) {
                for (int j=0; j<data.size(); j++) {
                    matrix[i][j] = data.get(j).getBins()[i];
                }
            }

            // now perform the sum
            final double[] output = new IndexRange(0, bins).mapToDouble(n -> MathUtils.sum(matrix[n]));

            // get the confidence interval, at the set width
            final ContaminationEstimate.ConfidenceInterval newInterval = new ContaminationEstimate.ConfidenceInterval(output);

            //TODO: ditto about a Table method
            logger.info(String.format("%s\t%s\t"+pm+"\t"+pm+"\t"+pm+"\t"+pm+"\t"+"%d",
                    aggregationKey,
                    "n/a",
                    newInterval.getContamination(),
                    (newInterval.getStop() - newInterval.getStart()),
                    newInterval.getStart(),
                    newInterval.getStop(),
                    data.size()));
        }
    }
}