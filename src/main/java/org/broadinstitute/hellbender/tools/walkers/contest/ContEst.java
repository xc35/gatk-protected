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
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.PrintStream;
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
 *     -I eval:tumor.bam \
 *     -I genotype:normal.bam \
 *     --popFile populationAlleleFrequencies.vcf \
 *     -L populationSites.interval_list
 *     [-L targets.interval_list] \
 *     -isr INTERSECTION \
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

    //TODO: what about Lane??
    public enum AggregationLevel {
        SAMPLE,    // calculate contamination for each sample
        READGROUP, // for each read group
        BAM        // for all inputs as a single source
    }

    // the population information; the allele frequencies for each position in known populations
    @Argument(fullName="popfile", shortName = "pf", doc="the variant file containing information about the population allele frequencies", optional = false)
    public FeatureInput<VariantContext> populationAlleleFrequencies;

    //FIXME - this is pretty clumsy because we don't have tagged inputs for bams yet
    @Argument(fullName = "evalSampleName", optional = false, doc = "eval sample name as in BAM header")
    public String evalSample;

    @Argument(fullName = "matchedNormalSampleName", optional = false, doc = "matched normal sample name as in BAM header")
    public String matchedNormalSample;

    // ------------------------------------------------------------------------------------------------------------------------------------------------------
    // outputs and args
    // ------------------------------------------------------------------------------------------------------------------------------------------------------

    @Argument(fullName = "min_qscore", optional = true, doc = "threshold for minimum base quality score")
    public int MIN_QSCORE = 20;

    @Argument(fullName = "min_mapq", optional = true, doc = "threshold for minimum mapping quality score")
    public int MIN_MAPQ = 20;

    @Argument(fullName = "beta_threshold", doc = "threshold for p(f>=0.5) to trim", optional = true)
    public double BETA_THRESHOLD = 0.95;

    @Argument(shortName = "agg", fullName = "aggregations", doc = "set to BAM (default), SAMPLE or READGROUP to produce per-bam, per-sample or per-lane estimates", optional = true)
    private Set<AggregationLevel> aggregations = null;

    @Argument(shortName = "pc", fullName = "precision", doc = "the degree of precision to which the contamination tool should estimate (e.g. the bin size)", optional = true)
    private double precision = 0.01;

    @Argument(shortName = "lf", fullName = "likelihood_file", doc = "write the likelihood values to the specified location", optional = true)
    public PrintStream likelihoodFile = null;

    @Argument(shortName = "population", fullName = "population", doc = "evaluate contamination for just a single contamination population", optional = true)
    public String population = "CEU";

    //TODO: populations should be enums
    private static final String[] ALL_POPULATIONS = {"ALL", "CHD", "LWK", "CHB", "CEU", "MXL", "GIH", "MKK", "TSI", "CLM", "GBR", "ASW", "YRI", "IBS", "FIN", "PUR", "JPT", "CHS"};

    // ------------------------------------------------------------------------------------------------------------------------------------------------------
    // global variables to the walker
    // ------------------------------------------------------------------------------------------------------------------------------------------------------
    private static final Allele[] ALLELES = {Allele.create((byte) 'A'), Allele.create((byte) 'C'), Allele.create((byte) 'G'), Allele.create((byte) 'T')};

    private final Map<String, AggregationLevel> contaminationNames = new LinkedHashMap<>();       // a list, containing the contamination names, be it read groups or bam file names

    Collection<String> allSamples;
    Collection<String> allReadGroups;

    //TODO:
    private String[] populationsToEvaluate;

    int countGenotypeHomVar = 0;

    ContaminationResults accumulatedResult = new ContaminationResults();


    @Override
    public void onTraversalStart() {
        verifySamples();

        if (aggregations == null) {
            aggregations = new LinkedHashSet<>();
            aggregations.add(AggregationLevel.BAM);
        }

        allReadGroups = getHeaderForReads().getReadGroups().stream().map(SAMReadGroupRecord::getId).collect(Collectors.toList());
        allSamples = getHeaderForReads().getReadGroups().stream().map(SAMReadGroupRecord::getSample).distinct().collect(Collectors.toList());

        this.populationsToEvaluate = (population == null || "EVERY".equals(population)) ? ALL_POPULATIONS : new String[]{population};

    }

    private void verifySamples() {
        final SampleList samplesList = new IndexedSampleList(new ArrayList<>(ReadUtils.getSamplesFromHeader(getHeaderForReads())));
        if (!samplesList.asListOfSamples().contains(evalSample)) {
            throw new UserException.BadInput("BAM header sample names " + samplesList.asListOfSamples() + "does not contain given tumor" +
                    " sample name " + evalSample);
        } else if (!samplesList.asListOfSamples().contains(matchedNormalSample)) {
            throw new UserException.BadInput("BAM header sample names " + samplesList.asListOfSamples() + "does not contain given normal" +
                    " sample name " + matchedNormalSample);
        }
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
        final List<VariantContext> vcs = features.getValues(populationAlleleFrequencies);
        if (vcs.isEmpty()) {
            return;
        }

        final VariantContext popVC = vcs.get(0);
        final Genotype genotype = getGenotypeFromMatchedNormal(context, ref, getHeaderForReads());

        // only use homozygous sites
        if (genotype == null || !genotype.isHomVar()) {
            return;
        } else {
            countGenotypeHomVar++;
        }

        // only use non-reference sites
        final byte myBase = genotype.getAllele(0).getBases()[0];

        // our map of contamination results
        final Map<String, Map<String, ContaminationStats>> contaminationResults = new LinkedHashMap<>();
        final ReadPileup wholePileup = context.getBasePileup().getPileupForSample(evalSample, getHeaderForReads());

        // if we're by-lane, get those stats
        for (final Map.Entry<String, AggregationLevel> namePair : contaminationNames.entrySet()) {
            final ReadPileup pile = getPileupForAggregationLevel(wholePileup, namePair);
            final ReadPileup filteredPile =
                    pile.makeFilteredPileup(pe -> pe.getQual() >= MIN_QSCORE && pe.getMappingQual() >= MIN_MAPQ);

            final Map<String, ContaminationStats> results = calcStats(
                            filteredPile,
                            myBase,
                            popVC,
                    populationsToEvaluate);

            if (!results.isEmpty()) {
                contaminationResults.put(namePair.getKey(), results);
            }
        }
        // add this site to our collected stats
        accumulatedResult.add(contaminationResults);
    }

    private ReadPileup getPileupForAggregationLevel(ReadPileup wholePileup, Map.Entry<String, AggregationLevel> namePair) {
        if (namePair.getValue() == AggregationLevel.READGROUP) {
            return wholePileup.makeFilteredPileup(pe -> Objects.equals(pe.getRead().getReadGroup(), namePair.getKey()));
        } else if (namePair.getValue() == AggregationLevel.SAMPLE) {
            return wholePileup.getPileupForSample(namePair.getKey(), getHeaderForReads());
        } else {
            return wholePileup;
        }
    }

    private Genotype getGenotypeFromMatchedNormal(final AlignmentContext context, final ReferenceContext referenceContext, final SAMFileHeader header) {
        final ReadPileup pileup = context.getBasePileup().getPileupForSample(matchedNormalSample, header);
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

    private static PopulationFrequencyInfo parsePopulationFrequencyInfo(final VariantContext variantContext, final String population) {
        PopulationFrequencyInfo info = null;

        @SuppressWarnings("unchecked")
        final List<String> values = (List<String>) variantContext.getAttribute(population);

        if (values != null) {
            byte majorAllele = 0;
            byte minorAllele = 0;
            double maf = -1;

            for (String str : values) {
                // strip off the curly braces and trim whitespace
                if (str.startsWith("{")) {
                    str = str.substring(1, str.length());
                }
                if (str.contains("}")) {
                    str = str.substring(0, str.indexOf("}"));
                }
                str = str.trim();
                final String[] spl = str.split("=");

                final byte allele = (byte) spl[0].trim().charAt(0);
                final double af = Double.valueOf(spl[1].trim());

                if (af <= 0.5 && minorAllele == 0) {
                    minorAllele = allele;
                    maf = af;
                } else {
                    majorAllele = allele;
                }

            }

            info = new PopulationFrequencyInfo(majorAllele, minorAllele, maf);
        }
        return info;
    }


    /**
     * Calculate the contamination values per division, be it lane, meta, sample, etc
     * @param pileup the pileup
     * @param myAllele the allele we have (our hom var genotype allele)
     * @param popVC the population variant context from hapmap
     * @param pops contaminating populations to run over
     * @return a mapping of each target population to their estimated contamination
     */
    //TODO: I hate that thsi operates via a side effect of printing
    private Map<String, ContaminationStats> calcStats(final ReadPileup pileup,
                                                      final byte myAllele,
                                                      final VariantContext popVC,
                                                      final String[] pops) {
        final Map<String, ContaminationStats> ret = new LinkedHashMap<>();

        for (final String pop : pops) {
            final PopulationFrequencyInfo info = parsePopulationFrequencyInfo(popVC, pop);
            final double alleleFreq = info.getMinorAlleleFrequency();
            if (alleleFreq > 0.5) {
                throw new RuntimeException("Minor allele frequency is greater than 0.5, this is an error; we saw AF of " + alleleFreq);
            }

            final long majorCounts = counter.get(Nucleotide.valueOf(info.getMajorAllele()));
            final long minorCounts = counter.get(Nucleotide.valueOf(info.getMinorAllele()));
            final long otherCounts = pileup.size() - majorCounts - minorCounts;

            // only use sites where this is the minor allele
            if (myAllele == info.minorAllele) {
                final ContaminationEstimate est = new ContaminationEstimate(precision, alleleFreq, pileup, info.getMinorAllele(), info.getMajorAllele(), pop);
                ret.put(pop, new ContaminationStats(1, alleleFreq, minorCounts, majorCounts, otherCounts, counter, est));
            }
        }
        return ret;
    }

    /**
     * Output all the stats to the appropriate files
     */
    @Override
    public Object onTraversalSuccess() {
        // filter out lanes / samples that don't have the minBaseCount
        final Map<String, Map<String, ContaminationStats>> cleanedMap = new LinkedHashMap<>();
        for (final Map.Entry<String, Map<String, ContaminationStats>> entry : accumulatedResult.getStats().entrySet()) {

            final Map<String, ContaminationStats> newMap = new LinkedHashMap<>();

            final Map<String, ContaminationStats> statMap = entry.getValue();
            for (final String popKey : statMap.keySet()) {
                final ContaminationStats stat = statMap.get(popKey);
                    newMap.put(popKey, stat);
            }

            cleanedMap.put(entry.getKey(), newMap);

        }

        // output results at the end, based on the input parameters
        accumulatedResult.setStats(cleanedMap);
        accumulatedResult.outputReport(precision, BETA_THRESHOLD);
        if (likelihoodFile != null) {
            accumulatedResult.writeCurves(likelihoodFile);
        }
        logger.info("Homozygous variant sites: " + countGenotypeHomVar);

        return accumulatedResult;
    }
}
