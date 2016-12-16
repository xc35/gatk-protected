package org.broadinstitute.hellbender.tools.exome.segmentation;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.primitives.Doubles;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.IntStream;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class CopyRatioSegmenter extends ScalarHMMSegmenter<Double> {
    private double logCoverageCauchyWidth;

    private static final double DEFAULT_INITIAL_CAUCHY_WIDTH = 0.1;
    private static final double NEUTRAL_LOG_2_COPY_RATIO = 0.0;
    private static final double MAX_REASONABLE_CAUCHY_WIDTH = 1.0;
    private static final double MIN_LOG_2_COPY_RATIO = -5.0;
    private static final double MAX_LOG_2_COPY_RATIO = 5.0;

    /**
     * @param initialNumStates  A liberal estimate of the number of hidden minor allele fraction values to
     *                          include in the model.  Hidden states are pruned as the model is learned.
     */
    public CopyRatioSegmenter(final int initialNumStates, final List<SimpleInterval> positions, final List<Double> data) {
        super(initialNumStates, positions, data, initialLog2CopyRatios(initialNumStates));
        logCoverageCauchyWidth = DEFAULT_INITIAL_CAUCHY_WIDTH;
    }

    /**
     * evenly-spaced log-2 copy ratios
     * @param K the initial number of hidden states
     */
    private static List<Double> initialLog2CopyRatios(final int K) {
        final List<Double> result = new ArrayList<>();
        result.add(NEUTRAL_LOG_2_COPY_RATIO);   // 2, neutral is never pruned and goes first
        result.add(MIN_LOG_2_COPY_RATIO);       // 0
        result.add(ParamUtils.log2(0.5));       // 1
        IntStream.range(3, K).forEach(n -> result.add(ParamUtils.log2(n / 2.0)));
        return result;
    }

    @Override
    protected void relearnHiddenStateValues(final ExpectationStep eStep) { }

    @Override
    protected ClusteringGenomicHMM<Double, Double> makeModel() {
        return new CopyRatioHiddenMarkovModel(Doubles.toArray(hiddenStateValues), weights, getMemoryLength(), logCoverageCauchyWidth);
    }

    @Override
    protected void relearnAdditionalParameters(final ExpectationStep eStep) {
        //relearn the Cauchy width of the emission distribution
        final Function<Double, Double> emissionLogLikelihood = width -> {
            double logLikelihood = 0.0;
            for (int position = 0; position < positions.size(); position++) {
                for (int state = 0; state < weights.length; state++) {
                    final double eStepPosterior = eStep.pStateAtPosition(state, position);
                    logLikelihood += eStepPosterior < NEGLIGIBLE_POSTERIOR_FOR_M_STEP ? 0 : eStepPosterior
                            * CopyRatioHiddenMarkovModel.logEmissionProbability(data.get(position), hiddenStateValues.get(state), width);
                }
            }
            return logLikelihood;
        };

        logCoverageCauchyWidth = OptimizationUtils.argmax(emissionLogLikelihood, 0, MAX_REASONABLE_CAUCHY_WIDTH, logCoverageCauchyWidth,
                RELATIVE_TOLERANCE_FOR_OPTIMIZATION, ABSOLUTE_TOLERANCE_FOR_OPTIMIZATION, MAX_EVALUATIONS_FOR_OPTIMIZATION);

        logger.info("New coverage standard deviation learned: " + logCoverageCauchyWidth);
    }

    @Override
    protected double minHiddenStateValue() { return MIN_LOG_2_COPY_RATIO; }

    @Override
    protected double maxHiddenStateValue() { return  MAX_LOG_2_COPY_RATIO; }

    public double getLogCoverageCauchyWidth() {
        return logCoverageCauchyWidth;
    }

}
