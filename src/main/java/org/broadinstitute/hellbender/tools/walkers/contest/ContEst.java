package org.broadinstitute.hellbender.tools.walkers.contest;

import htsjdk.variant.variantcontext.*;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.util.MathArrays;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

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
public final class ContEst extends VariantWalker {
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

    //pairs of AggregationLevel with their key eg "BAM" or the sample name or the read group ID
    final List<ImmutablePair<AggregationLevel, String>> aggregationKeys = new ArrayList<>();

    final MutableInt homVarCount = new MutableInt(0);

    // a map of our contamination targets and their cumulativeEstimates
    // key: aggregation entity ("BAM", sample name, or read group ID)
    // value: cumulative ContaminationEstimate
    private final Map<String, double[]> cumulativeLogLikelihoods = new HashMap<>();

    // the likelihood as a function of contamination is not analytic so we simply track it over a discrete array of values
    private double[] discreteContaminations;

    @Override
    public void onTraversalStart() {
        discreteContaminations = GATKProtectedMathUtils.createEvenlySpacedPoints(0, 1, (int) Math.round(1.0 / precision) + 1);
        final SampleList samplesList = new IndexedSampleList(new ArrayList<>(ReadUtils.getSamplesFromHeader(getHeaderForReads())));
        if (!samplesList.asListOfSamples().contains(evalSample)) {
            throw new UserException.BadInput("BAM header sample names " + samplesList.asListOfSamples() + "does not contain given tumor" + " sample name " + evalSample);
        }

        if (aggregations == null || aggregations.contains(AggregationLevel.BAM)) {
            aggregationKeys.add(new ImmutablePair<>(AggregationLevel.BAM, "BAM"));
        }
        if (aggregations.contains(AggregationLevel.SAMPLE)) {
            getHeaderForReads().getReadGroups().stream().map(rg -> new ImmutablePair<>(AggregationLevel.SAMPLE, rg.getSample())).distinct()
                    .forEach(aggregationKeys::add);
        }
        if (aggregations.contains(AggregationLevel.READGROUP)) {
            getHeaderForReads().getReadGroups().stream().map(rg -> new ImmutablePair<>(AggregationLevel.READGROUP, rg.getId()))
                    .forEach(aggregationKeys::add);
        }

        aggregationKeys.stream().map(ImmutablePair::getRight).forEach(
                key -> cumulativeLogLikelihoods.put(key, new double[discreteContaminations.length]));
    }

    /**
     * our map function, which emits a contamination stats for each of the subgroups (lanes, samples, etc) that we encounter
     *
     * @return a mapping of our subgroup name to contamination estimate
     */
    @Override
    public void apply(final VariantContext exacSite, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final int position = exacSite.getStart();
        final List<GATKRead> reads = StreamSupport.stream(readsContext.spliterator(), false).collect(Collectors.toList());

        //TODO: verify that vc.getStart() - read.getAssignedStart() is the correct offset for creating a Pileup
        final List<Integer> offsets = reads.stream().map(read -> position - read.getAssignedStart()).collect(Collectors.toList());
        final ReadPileup pileup = new ReadPileup(exacSite, reads, offsets)
                .makeFilteredPileup(pe -> pe.getQual() >= MIN_QSCORE && pe.getMappingQual() >= MIN_MAPQ);

        //TODO: what about pileup for just the eval sample?


        final Genotype genotype = getGenotype(pileup, referenceContext);

        // only use homozygous sites
        if (genotype == null || !genotype.isHomVar()) {
            return;
        } else {
            homVarCount.increment();
        }

        // only use non-reference sites
        final byte myBase = genotype.getAllele(0).getBases()[0];



        for (final ImmutablePair<AggregationLevel, String> aggregationKey : aggregationKeys) {
            final AggregationLevel level = aggregationKey.getLeft();
            final String key = aggregationKey.getRight();
            final ReadPileup pileupForAggregation = getPileupForAggregation(pileup, level, key);
            final double[] siteLogLikelihoods = calcStats(pileupForAggregation, myBase, exacSite);
            addSiteEstimate(key, siteLogLikelihoods);
        }
    }

    final ReadPileup getPileupForAggregation(final ReadPileup pileup, final AggregationLevel level, final String key) {
        if (level == AggregationLevel.BAM) {
            return pileup;  // the key is an arbitrary placeholder
        } else if (level == AggregationLevel.SAMPLE) {
            return pileup.getPileupForSample(key, getHeaderForReads()); // the key is a sample name
        } else if (level == AggregationLevel.READGROUP) {
            return pileup.makeFilteredPileup(pe -> pe.getRead().getReadGroup().equals(key));    // the key is a read group ID
        } else {
            throw new IllegalStateException("Programmer error: unrecognized aggregation level.");
        }
    }

    // siteLogLikelihoods can be null, in which case nothing happens
    private void addSiteEstimate(final String aggregationKey, final double[] siteLogLikelihoods) {
        if (siteLogLikelihoods != null) {
            double[] currentLikelihoods = cumulativeLogLikelihoods.get(aggregationKey);
            cumulativeLogLikelihoods.put(aggregationKey, MathArrays.ebeAdd(currentLikelihoods, siteLogLikelihoods));
        }
    }

    private Genotype getGenotype(final ReadPileup pileup, final ReferenceContext referenceContext) {
        //TODO: this is a really dumb way of calling genotypes.  Replace with a posterior probability.  Then it won't be
        //TODO: necessary to use a minimum depth, which I have already gotten rid of in anticipation
        final double MIN_GENOTYPE_RATIO = 0.8;
        final int[] bases = pileup.getBaseCounts();
        final int maxAlleleIndex = GATKProtectedMathUtils.maxIndex(bases);
        final byte majorBase = BaseUtils.baseIndexToSimpleBase(maxAlleleIndex);
        final int numReads = (int) MathUtils.sum(bases);
        final byte refBase = referenceContext.getBase();
        final boolean isHomAlt = majorBase != refBase && bases[maxAlleleIndex] / (double) numReads >= MIN_GENOTYPE_RATIO;
        return isHomAlt ? new GenotypeBuilder(evalSample, Collections.singletonList(Allele.create(majorBase))).make() : null;
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

    private static PopulationFrequencyInfo parseExACInfo(final VariantContext exacSite) {
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
    private double[] calcStats(final ReadPileup pileup, final byte myAllele, final VariantContext exacVC) {
            final PopulationFrequencyInfo info = parseExACInfo(exacVC);
            final double alleleFreq = info.getMinorAlleleFrequency();

            // only use sites where this is the minor allele
            if (myAllele == info.minorAllele) {
                return new ContaminationEstimate(precision, alleleFreq, pileup, info.getMinorAllele(), info.getMajorAllele());
            } else {
                return null;
            }
    }

    @Override
    public Object onTraversalSuccess() {
        for (final ImmutablePair<AggregationLevel, String> aggregationKey : aggregationKeys) {
            final AggregationLevel level = aggregationKey.getLeft();
            final String key = aggregationKey.getRight();
            final double[] logLikelihood = cumulativeLogLikelihoods.get(key);
            //TODO: get confidence interval

            //TODO: table output
        }
        logger.info("Homozygous variant sites: " + homVarCount);
        return "SUCCESS";
    }

    //TODO: move the following to an Engine class

    // log likelihood of base pileup given that our sample is homozygous for an uncontaminated allele and
    // is contaminated with a different base with a given population frequency
    // we assume that the uncontaminated allele's and contaminating allele's frequencies add up to 1 i.e. thhis site
    // is very nearly biallelic
    public static double contaminationLogLikelihood(final double contamination,
                                                    final double contaminantAlleleFrequency,
                                                    final byte[] bases,
                                                    final byte[] quals,
                                                    final byte sampleBase,
                                                    final byte contaminatingBase) {
        Utils.validateArg(0.0 <= contaminantAlleleFrequency && contaminantAlleleFrequency <= 1.0, () -> "Invalid allele Freq: must be between 0 and 1 (inclusive), maf was " + contaminantAlleleFrequency);
        Utils.validateArg(bases.length == quals.length, "Must have same number of bases and quals");

        // probability that the contaminating sample is homozygous for the contaminating allele, heterozygous
        // or heterozygous for the sample allele
        final double logPHomContaminating = Math.log(MathUtils.square(contaminantAlleleFrequency));
        final double logPHetContaminating = Math.log(2 * contaminantAlleleFrequency * (1 - contaminantAlleleFrequency));
        final double logPHomUncontaminating = Math.log(MathUtils.square(1 - contaminantAlleleFrequency));

        final DoubleUnaryOperator pileupLogLikelihood = genotypeFractionOfContaminatingAllele ->
                logLikelihood(contamination, bases, quals, sampleBase, contaminatingBase, genotypeFractionOfContaminatingAllele);

        return GATKProtectedMathUtils.logSumExp(logPHomContaminating + pileupLogLikelihood.applyAsDouble(1.0),
                logPHetContaminating + pileupLogLikelihood.applyAsDouble(0.5),
                logPHomUncontaminating + pileupLogLikelihood.applyAsDouble(0.0));
    }

    // the log likelihood of a pileup given the fraction -- 0, 1/2, or 1 -- of the contaminating allele in the
    // contaminant genotype
    private static double logLikelihood(double contamination, byte[] bases, byte[] quals, byte sampleBase, byte contaminatingBase, double genotypeFractionOfWrongAllele) {
        return new IndexRange(0, bases.length).sum(n -> logLikelihood(sampleBase, contaminatingBase,
                contamination * genotypeFractionOfWrongAllele, QualityUtils.qualToErrorProb(quals[n]), bases[n]));
    }

    // If sample is homozygous for sampleBase and a fraction contaminatingFraction of the library DNA at this locus
    // has contaminatingBase, return the log likelihood to sequence sequencedBase given an error probability
    private static double logLikelihood(byte sampleBase, byte contaminatingBase, double contaminatingFraction, double pError, byte sequencedBase) {
        if (sequencedBase == sampleBase) {
            // sample DNA without error + contaminant DNA with error
            return Math.log(( 1 - contaminatingFraction) * (1 - pError) + contaminatingFraction * pError/3);
        } else if (sequencedBase == contaminatingBase) {
            // sample DNA with error + contaminant DNA without error
            return Math.log((1 - contaminatingFraction) * pError/3 + contaminatingFraction * (1 - pError));
        } else {
            // another allele
            return 0.0;
        }
    }
}