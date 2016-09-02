package org.broadinstitute.hellbender.tools.exome.segmentation;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.annotations.VisibleForTesting;
import com.google.cloud.dataflow.sdk.repackaged.com.google.common.primitives.Doubles;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.DecileCollection;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by davidben on 10/5/16.
 */
@CommandLineProgramProperties(
        summary = "Segment genomic data into regions of constant copy ratio and allele fraction.  Only supports one sample input.",
        oneLineSummary = "Segment genomic data into regions of constant copy ratio and allele fraction",
        programGroup = CopyNumberProgramGroup.class
)
public class PerformJointSegmentation extends CommandLineProgram {
    protected static final String INITIAL_NUM_COPY_RATIO_STATES_LONG_NAME = "initialNumberOfCopyRatioStates";
    protected static final String INITIAL_NUM_COPY_RATIO_STATES_SHORT_NAME = "initialNumCRStates";

    protected static final String INITIAL_NUM_ALLELE_FRACTION_STATES_LONG_NAME = "initialNumberOfAlleleFractionStates";
    protected static final String INITIAL_NUM_ALLELE_FRACTION_STATES_SHORT_NAME = "initialNumAFStates";

    @Argument(
            doc = "Tangent-normalized log2 read counts file",
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            fullName =  ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            optional = false
    )
    protected String coverageFile;

    @Argument(
            doc = "Input file for tumor-sample ref/alt read counts at normal-sample heterozygous-SNP sites (output of GetHetCoverage tool).",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpCountsFile;

    @Argument(
            doc = "Initial number of hidden copy-ratio states",
            fullName = INITIAL_NUM_COPY_RATIO_STATES_LONG_NAME,
            shortName = INITIAL_NUM_COPY_RATIO_STATES_SHORT_NAME,
            optional = false
    )
    protected int initialNumCRStates;

    @Argument(
            doc = "Initial number of hidden allele-fraction states",
            fullName = INITIAL_NUM_ALLELE_FRACTION_STATES_LONG_NAME,
            shortName = INITIAL_NUM_ALLELE_FRACTION_STATES_SHORT_NAME,
            optional = false
    )
    protected int initialNumAFStates;

    @Argument(
            doc = "Input file for allelic-bias panel of normals.",
            fullName = ExomeStandardArgumentDefinitions.ALLELIC_PON_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.ALLELIC_PON_FILE_SHORT_NAME,
            optional = true
    )
    protected File allelicPoNFile;

    @Argument(
            doc = "Output file for copy-ratio segments.  Output is not logged.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            optional = false
    )
    protected File outputSegmentsFile;

    @Override
    public Object doWork() {
        final AllelicPanelOfNormals allelicPoN =
                allelicPoNFile != null ? AllelicPanelOfNormals.read(allelicPoNFile) : AllelicPanelOfNormals.EMPTY_PON;
        final AllelicCountCollection acc = new AllelicCountCollection(snpCountsFile);
        final AlleleFractionSegmenter alleleFractionSegmenter = new AlleleFractionSegmenter(initialNumAFStates, acc, allelicPoN);

        final ReadCountCollection counts;
        try {
            counts = ReadCountCollectionUtils.parse(new File(coverageFile));
        } catch (final IOException ex) {
            throw new UserException.BadInput("could not read input file");
        }

        final List<Double> coverage = Doubles.asList(counts.counts().getColumn(0));
        final List<SimpleInterval> targets = counts.targets().stream().map(Target::getInterval).collect(Collectors.toList());
        final CopyRatioSegmenter copyRatioSegmenter = new CopyRatioSegmenter(initialNumCRStates, targets, coverage);

        final JointAFCRSegmenter jointSegmenter = JointAFCRSegmenter.createJointSegmenter(copyRatioSegmenter, alleleFractionSegmenter);
        final List<Pair<SimpleInterval, AFCRHiddenState>> segmentation = jointSegmenter.findSegments();

        final List<ACNVModeledSegment> segments = segmentation.stream().map(pair ->
                new ACNVModeledSegment(pair.getLeft(), errorlessPosterior(pair.getRight().getLog2CopyRatio()), errorlessPosterior(pair.getRight().getMinorAlleleFraction())))
                .collect(Collectors.toList());

        SegmentUtils.writeACNVModeledSegmentFile(outputSegmentsFile, segments, new Genome(counts, acc.getCounts()));

        return "SUCCESS";
    }

    @VisibleForTesting
    protected static PosteriorSummary errorlessPosterior(final double value) {
        final PosteriorSummary result = new PosteriorSummary(value, value, value);
        result.setDeciles(new DecileCollection(Arrays.asList(value)));
        return result;
    }
}
