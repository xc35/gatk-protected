package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.exome.HashedListTargetCollection;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by davidben on 6/6/16.
 */
public final class CopyRatioSegmenterUnitTest {
    private static final Random RANDOM = new Random(271);
    @Test
    public void testSegmentation() {
        final double[] trueWeights = new double[] {0.2, 0.5, 0.3};
        final double[] trueLog2CopyRatios = new double[] {-2.0, 0.0, 1.7};
        final double trueMemoryLength = 1e5;
        final double trueStandardDeviation = 0.2;

        final CopyRatioHiddenMarkovModel trueModel = new CopyRatioHiddenMarkovModel(trueLog2CopyRatios, trueWeights,
                trueMemoryLength, trueStandardDeviation);

        final int chainLength = 10000;
        final List<SimpleInterval> positions = randomPositions("chr1", chainLength, RANDOM, trueMemoryLength/4);
        final List<Integer> trueStates = trueModel.generateHiddenStateChain(positions);
        final List<Double> trueCopyRatioSequence = trueStates.stream().map(n -> trueLog2CopyRatios[n]).collect(Collectors.toList());

        final List<Double> data = trueCopyRatioSequence.stream()
                .map(cr -> generateData(trueStandardDeviation, cr)).collect(Collectors.toList());

        final CopyRatioSegmenter segmenter = new CopyRatioSegmenter(10, positions, data);
        final List<Pair<SimpleInterval, Double>> segmentation = segmenter.findSegments();

        //TODO: extract this common motif
        final TargetCollection<SimpleInterval> tc = new HashedListTargetCollection<>(positions);
        final List<ModeledSegment> segments = segmentation.stream().map(pair ->
                new ModeledSegment(pair.getLeft(), tc.targetCount(pair.getLeft()), pair.getRight())).collect(Collectors.toList());

        final double[] segmentCopyRatios = segments.stream()
                .flatMap(s -> Collections.nCopies((int) s.getTargetCount(), s.getSegmentMean()).stream())
                .mapToDouble(x -> x).toArray();

        final double averageCopyRatioError = IntStream.range(0, trueCopyRatioSequence.size())
                .mapToDouble(n -> Math.abs(segmentCopyRatios[n] - trueCopyRatioSequence.get(n)))
                .average().getAsDouble();

        Assert.assertEquals(averageCopyRatioError, 0, 0.025);
    }

    protected static double generateData(double trueStandardDeviation, Double cr) {
        return cr + RANDOM.nextGaussian() * trueStandardDeviation;
    }

    @Test
    public void testChromosomesOnDifferentSegments() {
        final double[] trueLog2CopyRatios = new double[] {-2.0, 0.0, 1.7};
        final double trueMemoryLength = 1e5;

        final double trueStandardDeviation = 0.2;

        // randomly set positions
        final int chainLength = 100;
        final List<SimpleInterval> positions = randomPositions("chr1", chainLength, RANDOM, trueMemoryLength/4);
        positions.addAll(randomPositions("chr2", chainLength, RANDOM, trueMemoryLength/4));
        positions.addAll(randomPositions("chr3", chainLength, RANDOM, trueMemoryLength/4));

        final int trueState = 2;    //fix everything to the same state 2

        final List<Double> data = new ArrayList<>();
        for (int n = 0; n < positions.size(); n++) {
            final double copyRatio = trueLog2CopyRatios[trueState];
            final double observed = generateData(trueStandardDeviation, copyRatio);
            data.add(observed);
        }

        final CopyRatioSegmenter segmenter = new CopyRatioSegmenter(10, positions, data);
        final List<Pair<SimpleInterval, Double>> segmentation = segmenter.findSegments();
        final TargetCollection<SimpleInterval> tc = new HashedListTargetCollection<>(positions);
        final List<ModeledSegment> segments = segmentation.stream().map(pair ->
                new ModeledSegment(pair.getLeft(), tc.targetCount(pair.getLeft()), Math.pow(2, pair.getRight()))).collect(Collectors.toList());

        //check that each chromosome has at least one segment
        final int numDifferentContigsInSegments = (int) segments.stream().map(ModeledSegment::getContig).distinct().count();
        Assert.assertEquals(numDifferentContigsInSegments, 3);
    }

    public static List<SimpleInterval> randomPositions(final String contig, final int chainLength, final Random random, final double separationScale) {
        final List<SimpleInterval> positions = new ArrayList<>();
        int position = 1;
        for (int n = 0; n < chainLength; n++) {
            position += random.nextInt((int) separationScale);
            final SimpleInterval interval = new SimpleInterval(contig, position, position);
            positions.add(interval);
        }
        return positions;
    }
}