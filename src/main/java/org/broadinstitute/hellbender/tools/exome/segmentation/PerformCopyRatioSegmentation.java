package org.broadinstitute.hellbender.tools.exome.segmentation;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.primitives.Doubles;
import org.apache.avro.test.Simple;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by davidben on 5/23/16.
 */
@CommandLineProgramProperties(
        summary = "Segment genomic data into regions of constant copy ratio.  Only supports one sample input.",
        oneLineSummary = "Segment genomic data into regions of constant copy ratio",
        programGroup = CopyNumberProgramGroup.class
)
public final class PerformCopyRatioSegmentation extends CommandLineProgram {
    protected static final String INITIAL_NUM_STATES_LONG_NAME = "initialNumberOfStates";
    protected static final String INITIAL_NUM_STATES_SHORT_NAME = "initialNumStates";

    @Argument(
            doc = "Tangent-normalized log2 read counts file",
            fullName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected String coverageFile;

    @Argument(
            doc = "Initial number of hidden copy-ratio states",
            fullName = INITIAL_NUM_STATES_LONG_NAME,
            shortName = INITIAL_NUM_STATES_SHORT_NAME,
            optional = false
    )
    protected int initialNumStates;

    @Argument(
            doc = "Output file for copy-ratio segments.  Output is not logged.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            optional = false
    )
    protected File outputSegmentsFile;

    @Override
    public Object doWork() {
        final String sampleName = ReadCountCollectionUtils.getSampleNameForCLIsFromReadCountsFile(new File(coverageFile));
        final ReadCountCollection counts;
        try {
            counts = ReadCountCollectionUtils.parse(new File(coverageFile));
        } catch (final IOException ex) {
            throw new UserException.BadInput("could not read input file");
        }

        final List<Double> coverage = Doubles.asList(counts.counts().getColumn(0));
        final List<SimpleInterval> intervals = counts.targets().stream().map(Target::getInterval).collect(Collectors.toList());
        final CopyRatioSegmenter segmenter = new CopyRatioSegmenter(initialNumStates, intervals, coverage);
        final List<Pair<SimpleInterval, Double>> intervalLog2CoveragePairs = segmenter.findSegments();
        final TargetCollection<SimpleInterval> tc = new HashedListTargetCollection<>(intervals);

        final List<ModeledSegment> segments = intervalLog2CoveragePairs.stream().map(pair ->
                new ModeledSegment(pair.getLeft(), tc.targetCount(pair.getLeft()), Math.pow(2, pair.getRight()))).collect(Collectors.toList());
        SegmentUtils.writeModeledSegmentFile(outputSegmentsFile, segments, sampleName, true);

        return "SUCCESS";
    }
}
