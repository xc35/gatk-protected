package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.io.File;
import java.util.*;

@CommandLineProgramProperties(
        summary = "Filter somatic SNPs and indels called by Mutect2",
        oneLineSummary = "Filter somatic SNPs and indels called by Mutect2",
        programGroup = VariantProgramGroup.class
)
public final class FilterMutectCalls extends VariantWalker {

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output filtered VCF file", optional=false)
    private final String outputVcf = null;

    @ArgumentCollection
    protected M2FiltersArgumentCollection MTFAC = new M2FiltersArgumentCollection();

    private VariantContextWriter vcfWriter;

    private final List<VariantContext> resourceFilteredVcs = new ArrayList<>();

    private Mutect2FilteringEngine filteringEngine;


    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> headerLines = new HashSet<>(inputHeader.getMetaDataInSortedOrder());
        Mutect2FilteringEngine.M_2_FILTER_NAMES.stream().map(GATKVCFHeaderLines::getFilterLine).forEach(headerLines::add);
        headerLines.add(new VCFFilterHeaderLine(Mutect2FilteringEngine.ARTIFACT_IN_NORMAL_FILTER_NAME, "artifact_in_normal"));
        headerLines.add(new VCFFilterHeaderLine(Mutect2FilteringEngine.CONTAMINATION_FILTER_NAME, "contamination"));

        final VCFHeader vcfHeader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
        vcfWriter = createVCFWriter(new File(outputVcf));
        vcfWriter.writeHeader(vcfHeader);

        final String tumorSample = getHeaderForVariants().getMetaDataLine(Mutect2Engine.TUMOR_SAMPLE_KEY_IN_VCF_HEADER).getValue();
        filteringEngine = new Mutect2FilteringEngine(MTFAC, tumorSample, getBestAvailableSequenceDictionary());
    }

    @Override
    public Object onTraversalSuccess() {
        // TODO: implement sophisticated filtering
        for (final VariantContext vc : resourceFilteredVcs) {
            final VariantContextBuilder vcb = new VariantContextBuilder(vc);
            vcb.filters(filteringEngine.calculateAnnotationBasedFilters(vc));
            vcfWriter.add(vcb.make());
        }
        return "SUCCESS";
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {
        final Set<String> resourceBasedFilters = filteringEngine.calculateResourceBasedFilters(vc, fc);
        final VariantContext resourceFilteredVc = resourceBasedFilters.isEmpty() ? vc : new VariantContextBuilder(vc).filters(resourceBasedFilters).make();
        resourceFilteredVcs.add(resourceFilteredVc);
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
