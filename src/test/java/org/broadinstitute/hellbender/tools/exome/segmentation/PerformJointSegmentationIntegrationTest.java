package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

import static org.testng.Assert.*;

/**
 * Created by davidben on 10/21/16.
 */
public class PerformJointSegmentationIntegrationTest extends CommandLineProgramTest {
    private static final String TOOLS_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";
    private static final File ALLELIC_COUNTS_FILE = new File(TOOLS_TEST_DIRECTORY, "snps-for-allelic-integration.tsv");
    private static final File LOG2_TN_COVERAGE_FILE = new File(TOOLS_TEST_DIRECTORY, "coverages-for-copy-ratio-modeller.tsv" );

    @Test
    public void testCommandLine() throws IOException {
        final File tnCoverageFile = LOG2_TN_COVERAGE_FILE;
        final File snpFile = ALLELIC_COUNTS_FILE;
        final File outputSegmentFile = createTempFile("segments", ".seg");
        final int initialNumCRStates = 10;
        final int initialNumAFStates = 10;
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, tnCoverageFile.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME, snpFile.getAbsolutePath(),
                "-" + PerformJointSegmentation.INITIAL_NUM_COPY_RATIO_STATES_SHORT_NAME, Integer.toString(initialNumCRStates),
                "-" + PerformJointSegmentation.INITIAL_NUM_ALLELE_FRACTION_STATES_SHORT_NAME, Integer.toString(initialNumAFStates),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, outputSegmentFile.getAbsolutePath(),
        };
        runCommandLine(arguments);

        final List<ACNVModeledSegment> segments = SegmentUtils.readACNVModeledSegmentFile(outputSegmentFile);
    }
}