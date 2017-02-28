package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredIntervalArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationRecord;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class M2FiltersArgumentCollection extends AssemblyBasedCallerArgumentCollection {
    private static final long serialVersionUID = 9345L;

    /***************************************/
    // Reference Metadata inputs
    /***************************************/
    /**
     * Mutect2 has the ability to use COSMIC data in conjunction with dbSNP to adjust the threshold for evidence of a variant
     * in the normal.  If a variant is present in dbSNP, but not in COSMIC, then more evidence is required from the normal
     * sample to prove the variant is not present in germline.
     */
    @Argument(fullName="cosmic", shortName = "cosmic", doc="VCF of COSMIC sites", optional = true)
    public FeatureInput<VariantContext> cosmic;

    /**
     * a vcf, for example a sites vcf from ExAC, containing exonic germline variants and their allele frequencies
     */
    @Argument(fullName="exomeResource", shortName = "exome", doc="VCF of exonic variants", optional = true)
    public FeatureInput<VariantContext> exomeResource;

    /**
     * a vcf, for example a sites vcf from ExAC, containing exonic germline variants and their allele frequencies
     */
    @Argument(fullName="alleleFrequencyOfVariantsNotInExomeResource", doc="Assumed allele frequency of exonic variants not found in exme resource", optional = true)
    public double afOfVariantsNotInExomeResource = 1.0 / 200_000;    //this is 1 / (size of ExAC)

    /**
     * Intervals of exome resource.  Useful because if a variant within these intervals is *not* in the exome resource
     * it tells us that the prior probability of a germline variant is very small.
     */
    @ArgumentCollection
    public IntervalArgumentCollection exomeIntervals = new IntervalArgumentCollection() {

        @Argument(fullName = "exomeIntervals", doc = "exome intervals corresponding to exome resource", optional = true)
        protected final List<String> intervalStrings = new ArrayList<>();

        @Override
        protected List<String> getIntervalStrings() {
            return intervalStrings;
        }

        @Override
        protected void addToIntervalStrings(String newInterval) {
            if ( traversalParameters != null ) {
                throw new IllegalStateException("addToIntervalStrings() cannot be called after interval parsing is complete");
            }

            intervalStrings.add(newInterval);
        }
    };


    /**
     * a vcf, for example a sites vcf from 1000 Genomes, containing whole-genome germline variants and their allele frequencies
     */
    @Argument(fullName="genome_resource", shortName = "genome", doc="VCF of whole-genome variants", optional = true)
    public FeatureInput<VariantContext> genomeResource;

    /**
     * a vcf, for example a sites vcf from ExAC, containing exonic germline variants and their allele frequencies
     */
    @Argument(fullName="alleleFrequencyOfVariantsNotInGenomeResource", doc="Assumed allele frequency of variants not found in wgs resource", optional = true)
    public double afOfVariantsNotInGenomeResource = 1.0 / 3_000;    //this is 1 / (size of 1000 Genomes)

    /**
     * A panel of normals can be a useful (optional) input to help filter out commonly seen sequencing noise that may appear as low allele-fraction somatic variants.
     */
    @Argument(fullName="normal_panel", shortName = "PON", doc="VCF file of sites observed in normal", optional = true)
    public FeatureInput<VariantContext> panelOfNormals;

    /**
     * This is a measure of the minimum evidence to support that a variant observed in the tumor is not also present in the normal
     * as an artifact i.e. not as a germline event.
     */
    @Argument(fullName = "normal_artifact_lod", optional = true, doc = "LOD threshold for calling normal artifacts")
    public double NORMAL_ARTIFACT_LOD_THRESHOLD = 0.0;

    /**
     * The LOD threshold for the normal is typically made more strict if the variant has been seen in dbSNP (i.e. another
     * normal sample). We thus require MORE evidence that a variant is NOT seen in this tumor's normal if it has been observed as a germline variant before.
     */
    @Argument(fullName = "dbsnpAlleleFrequency", optional = true, doc = "Assumed allele frequency of variants found in dbSNP")
    public double dbsnpAlleleFrequency = 1.0 / 100;

    @Hidden
    @Argument(fullName = "strand_artifact_lod", optional = true, doc = "LOD threshold for calling strand bias")
    public float STRAND_ARTIFACT_LOD_THRESHOLD = 2.0f;

    @Hidden
    @Argument(fullName = "strand_artifact_power_threshold", optional = true, doc = "power threshold for calling strand bias")
    public float STRAND_ARTIFACT_POWER_THRESHOLD = 0.9f;

    @Argument(fullName = "enable_strand_artifact_filter", optional = true, doc = "turn on strand artifact filter")
    public boolean ENABLE_STRAND_ARTIFACT_FILTER = true;

    @Argument(fullName = "enable_clustered_read_position_filter", optional = true, doc = "turn on clustered read position filter")
    public boolean ENABLE_CLUSTERED_READ_POSITION_FILTER = true;

    /**
     * This argument is used for the M1-style read position filter
     */
    @Argument(fullName = "pir_median_threshold", optional = true, doc="threshold for clustered read position artifact median")
    public double PIR_MEDIAN_THRESHOLD = 10;

    /**
     * This argument is used for the M1-style read position filter
     */
    @Argument(fullName = "pir_mad_threshold", optional = true, doc="threshold for clustered read position artifact MAD")
    public double PIR_MAD_THRESHOLD = 3;

    @Argument(shortName = "contaminationTable", fullName = "contaminationTable", optional = true, doc="Table containing contamination information.")
    public File contaminationTable = null;

}
