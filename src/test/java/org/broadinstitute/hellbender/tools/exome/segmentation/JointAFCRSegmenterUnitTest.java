package org.broadinstitute.hellbender.tools.exome.segmentation;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.primitives.Doubles;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.tools.exome.HashedListTargetCollection;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionGlobalParameters;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.testng.Assert.*;

/**
 * Created by davidben on 10/6/16.
 */
public class JointAFCRSegmenterUnitTest {

    @Test
    public void testSegmentation() {
        final double hetProportion = 0.25; // probability that a datum is a het i.e. #hets / (#hets + #targets)
        final double[] trueWeights = new double[] {0.2, 0.5, 0.3};
        final double[] trueMinorAlleleFractions = new double[] {0.12, 0.32, 0.5};
        final double[] trueLog2CopyRatios = new double[] {-2.0, 0.0, 1.7};
        final List<AFCRHiddenState> trueJointStates = IntStream.range(0, trueLog2CopyRatios.length)
                .mapToObj(n -> new AFCRHiddenState(trueMinorAlleleFractions[n], trueLog2CopyRatios[n]))
                .collect(Collectors.toList());
        final double trueMemoryLength = 1e5;
        final double trueCauchyWidth = 0.2;
        final AlleleFractionGlobalParameters trueAFParams = new AlleleFractionGlobalParameters(1.0, 0.01, 0.01);


        final JointAFCRHiddenMarkovModel trueJointModel = new JointAFCRHiddenMarkovModel(trueJointStates, trueWeights, trueMemoryLength,
                trueAFParams, AllelicPanelOfNormals.EMPTY_PON, trueCauchyWidth);

        // generate joint truth
        final Random random = new Random(271);
        final int chainLength = 10000;
        final List<SimpleInterval> positions = CopyRatioSegmenterUnitTest.randomPositions("chr1", chainLength, random, trueMemoryLength/4);
        final List<Integer> trueHiddenStates = trueJointModel.generateHiddenStateChain(positions);
        final List<AFCRHiddenState> trueAFCRSequence = trueHiddenStates.stream().map(trueJointModel::getHiddenStateValue).collect(Collectors.toList());
        final double[] trueCopyRatioSequence = trueAFCRSequence.stream().mapToDouble(AFCRHiddenState::getLog2CopyRatio).toArray();
        final double[] trueAlleleFractionSequence = trueAFCRSequence.stream().mapToDouble(AFCRHiddenState::getMinorAlleleFraction).toArray();

        // generate separate af and cr data
        final GammaDistribution biasGenerator = AlleleFractionSegmenterUnitTest.getGammaDistribution(trueAFParams);
        final double outlierProbability = trueAFParams.getOutlierProbability();
        final AllelicCountCollection afData = new AllelicCountCollection();
        final List<Double> crData = new ArrayList<>();
        final List<SimpleInterval> crPositions = new ArrayList<>();
        for (int n = 0; n < positions.size(); n++) {
            final SimpleInterval position = positions.get(n);
            final AFCRHiddenState jointState = trueAFCRSequence.get(n);
            final double minorFraction = jointState.getMinorAlleleFraction();
            final double log2CopyRatio = jointState.getLog2CopyRatio();

            if (random.nextDouble() < hetProportion) {  // het datum
                afData.add(AlleleFractionSegmenterUnitTest.generateAllelicCount(minorFraction, position, random, biasGenerator, outlierProbability));
            } else {    //target datum
                crPositions.add(position);
                crData.add(CopyRatioSegmenterUnitTest.generateData(trueCauchyWidth, log2CopyRatio));
            }
        }

        // make the joint segmenter
        final CopyRatioSegmenter crSegmenter = new CopyRatioSegmenter(20, crPositions, crData);
        final AlleleFractionSegmenter afSegmenter = new AlleleFractionSegmenter(20, afData, AllelicPanelOfNormals.EMPTY_PON);
        final JointAFCRSegmenter segmenter = JointAFCRSegmenter.createJointSegmenter(crSegmenter, afSegmenter);



        final TargetCollection<SimpleInterval> tc = new HashedListTargetCollection<>(positions);
        final List<Pair<SimpleInterval, AFCRHiddenState>> segmentation = segmenter.findSegments();
        final List<ACNVModeledSegment> jointSegments = segmentation.stream()
                .map(pair -> {
                    final SimpleInterval position = pair.getLeft();
                    final AFCRHiddenState jointState = pair.getRight();
                    final PosteriorSummary crSummary = PerformJointSegmentation.errorlessPosterior(jointState.getLog2CopyRatio());
                    final PosteriorSummary afSummary = PerformJointSegmentation.errorlessPosterior(jointState.getMinorAlleleFraction());
                    return new ACNVModeledSegment(position, crSummary, afSummary);
                })
                .collect(Collectors.toList());

        final double[] segmentCopyRatios = jointSegments.stream()
                .flatMap(s -> Collections.nCopies(tc.targetCount(s.getInterval()), s.getSegmentMeanPosteriorSummary().getCenter()).stream())
                .mapToDouble(x->x).toArray();
        final double[] segmentMinorFractions = jointSegments.stream()
                .flatMap(s -> Collections.nCopies(tc.targetCount(s.getInterval()), s.getMinorAlleleFractionPosteriorSummary().getCenter()).stream())
                .mapToDouble(x->x).toArray();

        final double averageMinorFractionError = Arrays.stream(MathArrays.ebeSubtract(trueAlleleFractionSequence, segmentMinorFractions))
                .map(Math::abs).average().getAsDouble();
        final double averageCopyRatioError = Arrays.stream(MathArrays.ebeSubtract(trueCopyRatioSequence, segmentCopyRatios))
                .map(Math::abs).average().getAsDouble();

        Assert.assertEquals(averageMinorFractionError, 0, 0.03);
        Assert.assertEquals(averageCopyRatioError, 0, 0.03);
        int j = 6;
    }

}