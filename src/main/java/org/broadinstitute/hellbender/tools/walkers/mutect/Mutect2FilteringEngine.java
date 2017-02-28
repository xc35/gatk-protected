package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.tools.exome.HashedListTargetCollection;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationRecord;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * Created by David Benjamin on 9/15/16.
 */
public class Mutect2FilteringEngine {

    public final static String ARTIFACT_IN_NORMAL_FILTER_NAME = "artifact_in_normal";
    public final static String CONTAMINATION_FILTER_NAME = "contamination";

    public static final List<String> M_2_FILTER_NAMES = Arrays.asList(GATKVCFConstants.STR_CONTRACTION_FILTER_NAME, GATKVCFConstants.PON_FILTER_NAME,
            GATKVCFConstants.HOMOLOGOUS_MAPPING_EVENT_FILTER_NAME, GATKVCFConstants.CLUSTERED_EVENTS_FILTER_NAME,
            GATKVCFConstants.TUMOR_LOD_FILTER_NAME, GATKVCFConstants.GERMLINE_RISK_FILTER_NAME, GATKVCFConstants.TRIALLELIC_SITE_FILTER_NAME,
            GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME);

    private final M2FiltersArgumentCollection MTFAC;
    private final double contamination;
    private final String tumorSample;
    private final List<Pair<Predicate<VariantContext>, String>> filters;
    private final TargetCollection<SimpleInterval> exomeIntervals;

    public Mutect2FilteringEngine(final M2FiltersArgumentCollection MTFAC, final String tumorSample, final SAMSequenceDictionary sequenceDict) {
        this.MTFAC = MTFAC;
        contamination = MTFAC.contaminationTable == null ? 0.0 : ContaminationRecord.readContaminationTable(MTFAC.contaminationTable).get(0).getContamination();
        this.tumorSample = tumorSample;

        filters = Arrays.asList(
                new ImmutablePair<>(vc -> failsEventDistanceFilter(vc), GATKVCFConstants.HOMOLOGOUS_MAPPING_EVENT_FILTER_NAME),
                new ImmutablePair<>(vc -> failsTriallelicFilter(vc), GATKVCFConstants.TRIALLELIC_SITE_FILTER_NAME),
                new ImmutablePair<>(vc -> failsArtifactInNormalFilter(vc), ARTIFACT_IN_NORMAL_FILTER_NAME),
                new ImmutablePair<>(vc -> failsClusteredReadPositionFilter(vc), GATKVCFConstants.CLUSTERED_EVENTS_FILTER_NAME),
                new ImmutablePair<>(vc -> failsArtifactInNormalFilter(vc), ARTIFACT_IN_NORMAL_FILTER_NAME),
                new ImmutablePair<>(vc -> failsStrandBiasFilter(vc), GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME),
                new ImmutablePair<>(vc -> failsSTRFilter(vc), GATKVCFConstants.STR_CONTRACTION_FILTER_NAME),
                new ImmutablePair<>(vc -> failsContaminationFilter(vc), CONTAMINATION_FILTER_NAME));

        exomeIntervals = MTFAC.exomeIntervals.intervalsSpecified() ?
                new HashedListTargetCollection<>(MTFAC.exomeIntervals.getIntervals(sequenceDict)) : null;
    }

    // very naive M1-style contamination filter -- remove calls with AF less than the contamination fraction
    private boolean failsContaminationFilter(final VariantContext vc) {
        final Genotype tumorGenotype = vc.getGenotype(tumorSample);
        final double[] alleleFractions = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(tumorGenotype, VCFConstants.ALLELE_FREQUENCY_KEY,
                () -> new double[] {1.0}, 1.0);
        return MathUtils.arrayMax(alleleFractions) < contamination;
    }

    private static boolean failsTriallelicFilter(final VariantContext vc) {
        return vc.getNAlleles() > 2;
    }

    private static boolean failsSTRFilter(final VariantContext vc) {
        // STR contractions, such as ACTACTACT -> ACTACT, are overwhelmingly false positives so we hard filter by default
        if (vc.isIndel()) {
            final int[] rpa = vc.getAttributeAsList(GATKVCFConstants.REPEATS_PER_ALLELE_KEY).stream()
                    .mapToInt(o -> Integer.parseInt(String.valueOf(o))).toArray();
            final String ru = vc.getAttributeAsString(GATKVCFConstants.REPEAT_UNIT_KEY, "");
            if (rpa != null && rpa.length > 1 && ru.length() > 1) {
                if (rpa[0] - rpa[1] == 1) {
                    return true;
                }
            }
        }
        return false;
    }

    private boolean failsClusteredReadPositionFilter(final VariantContext vc) {
        if (!MTFAC.ENABLE_CLUSTERED_READ_POSITION_FILTER) {
            return false;
        }
        final Double tumorFwdPosMedian = (Double) vc.getAttribute(GATKVCFConstants.MEDIAN_LEFT_OFFSET_KEY);
        final Double tumorRevPosMedian = (Double) vc.getAttribute(GATKVCFConstants.MEDIAN_RIGHT_OFFSET_KEY);
        final Double tumorFwdPosMAD = (Double) vc.getAttribute(GATKVCFConstants.MAD_MEDIAN_LEFT_OFFSET_KEY);
        final Double tumorRevPosMAD = (Double) vc.getAttribute(GATKVCFConstants.MAD_MEDIAN_RIGHT_OFFSET_KEY);
        //If the variant is near the read end (median threshold) and the positions are very similar (MAD threshold) then filter
        return (tumorFwdPosMedian != null && tumorFwdPosMedian <= MTFAC.PIR_MEDIAN_THRESHOLD
                && tumorFwdPosMAD != null && tumorFwdPosMAD <= MTFAC.PIR_MAD_THRESHOLD) ||
                (tumorRevPosMedian != null && tumorRevPosMedian <= MTFAC.PIR_MEDIAN_THRESHOLD
                        && tumorRevPosMAD != null && tumorRevPosMAD <= MTFAC.PIR_MAD_THRESHOLD);

    }

    // filter out anything called in tumor that would also be called in the normal if it were treated as a tumor.
    // this handles shared artifacts, such as ones due to alignment and any shared aspects of sequencing
    private boolean failsArtifactInNormalFilter(final VariantContext vc) {
        if (!( vc.hasAttribute(SomaticGenotypingEngine.NORMAL_ARTIFACT_LOD_ATTRIBUTE)
                && vc.hasAttribute(GATKVCFConstants.TUMOR_LOD_KEY))) {
            return false;
        }

        final double[] normalArtifactLods = getArrayAttribute(vc, SomaticGenotypingEngine.NORMAL_ARTIFACT_LOD_ATTRIBUTE);
        final double[] tumorLods = getArrayAttribute(vc, GATKVCFConstants.TUMOR_LOD_KEY);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);

        return normalArtifactLods[indexOfMaxTumorLod] > MTFAC.NORMAL_ARTIFACT_LOD_THRESHOLD;
    }

    private boolean failsStrandBiasFilter(final VariantContext vc) {
        if ( !MTFAC.ENABLE_STRAND_ARTIFACT_FILTER) {
            return false;
        }
        // set defaults to high LOD, zero power, which passes the filter
        final double forwardLod = vc.getAttributeAsDouble(GATKVCFConstants.TLOD_FWD_KEY, Double.POSITIVE_INFINITY);
        final double reverseLod = vc.getAttributeAsDouble(GATKVCFConstants.TLOD_REV_KEY, Double.POSITIVE_INFINITY);
        final double forwardPower = vc.getAttributeAsDouble(GATKVCFConstants.TUMOR_SB_POWER_FWD_KEY, 0.0);
        final double reversePower = vc.getAttributeAsDouble(GATKVCFConstants.TUMOR_SB_POWER_REV_KEY, 0.0);
        return (forwardPower > MTFAC.STRAND_ARTIFACT_POWER_THRESHOLD && forwardLod < MTFAC.STRAND_ARTIFACT_LOD_THRESHOLD) ||
                (reversePower > MTFAC.STRAND_ARTIFACT_POWER_THRESHOLD && reverseLod < MTFAC.STRAND_ARTIFACT_LOD_THRESHOLD);
    }


    private static boolean failsEventDistanceFilter(final VariantContext vc) {
        final Integer eventCount = vc.getAttributeAsInt(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY, -1);
        return eventCount > 2;
    }

    public Set<String> calculateAnnotationBasedFilters(final VariantContext vc) {
        return filters.stream()
                .filter(p -> p.getLeft().test(vc))
                .map(p -> p.getRight())
                .collect(Collectors.toSet());
    }

    private double calculateGermlinePriorProbability(final VariantContext vc, final FeatureContext fc) {
            final List<VariantContext> exomeVariants = fc.getValues(MTFAC.exomeResource);
            final List<VariantContext> genomeVariants = fc.getValues(MTFAC.genomeResource);

            /**
             * First handle the cases when the variant *is* in one of our resources.  We give the exome resource highest
             * priority (since it's almost certainly the biggest with the best-resolved allele frequency), followed by
             * the genome
             */
            if (!exomeVariants.isEmpty()) { //TODO: replace isEmpty with check for allele match
                final VariantContext exomeVc = exomeVariants.get(0);
                //TODO: magic constant default -- shouldn't matter since resource shoul dhave AF, but still
                return exomeVc.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, 0.001);
            } else if ( !genomeVariants.isEmpty() ){ //TODO: replace isEmpty with check for allele match
                final VariantContext genomeVc = genomeVariants.get(0);
                //TODO: magic constant default -- shouldn't matter since resource shoul dhave AF, but still
                return genomeVc.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, 0.001);
            }

        /**
         * Now handle the case of when the variant is missing from one of our resources.  First we handle the case of
         * an exonic variant that's not in our exome resource.  For non-exonic variants missing from the genome resource
         * we further check whether the variant appears in dbSNP
         */
            if (MTFAC.exomeResource != null && MTFAC.exomeIntervals.intervalsSpecified() && exomeIntervals.targetCount(vc) > 0) {
                return MTFAC.afOfVariantsNotInExomeResource;
            } else {    //TODO: replace isEmpty with check for allele match
                return fc.getValues(MTFAC.dbsnp.dbsnp).isEmpty() ? MTFAC.afOfVariantsNotInGenomeResource : MTFAC.dbsnpAlleleFrequency;
            }

    }

    // filters based on a panel of normals and vcfs of known germline variants -- these are implemented as look-ups
    // that can be streamed
    public Set<String> calculateResourceBasedFilters(final VariantContext vc, final FeatureContext fc) {






        final Set<String> filters = new HashSet<>();

        final boolean inPanelOfNormals = !fc.getValues(MTFAC.panelOfNormals).isEmpty();
        if (inPanelOfNormals) {
            filters.add(GATKVCFConstants.PON_FILTER_NAME);
        }

        final boolean inCosmic = !fc.getValues(MTFAC.cosmic).isEmpty();
        final boolean inDbsnp = !fc.getValues(MTFAC.dbsnp.dbsnp).isEmpty();
        final List<VariantContext> exomeVariants = fc.getValues(MTFAC.exomeResource);
        final List<VariantContext> genomeVariants = fc.getValues(MTFAC.genomeResource);

        return filters;


            final double[] tumorLods = getArrayAttribute(vc, GATKVCFConstants.TUMOR_LOD_KEY);
            final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);

            if (inDbsnp && !inCosmic ) {
                // take the normal LOD of the best somatic alt allele
                final double normalLod = getArrayAttribute(vc, GATKVCFConstants.NORMAL_LOD_KEY)[indexOfMaxTumorLod];
                if (normalLod < MTFAC.NORMAL_DBSNP_LOD_THRESHOLD) {
                    filters.add(GATKVCFConstants.GERMLINE_RISK_FILTER_NAME);
                }
            }
    }

    private static double[] getArrayAttribute(final VariantContext vc, final String attribute) {
        return GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, attribute, () -> null, -1);
    }


}
