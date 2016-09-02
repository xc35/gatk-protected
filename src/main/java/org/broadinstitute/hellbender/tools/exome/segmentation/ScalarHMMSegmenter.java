package org.broadinstitute.hellbender.tools.exome.segmentation;

import com.google.cloud.dataflow.sdk.options.Hidden;
import com.google.cloud.dataflow.sdk.repackaged.com.google.common.primitives.Doubles;
import org.broadinstitute.hellbender.utils.*;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * HMM segmenter for scalar hidden data, such as allele fraction and copy ratio, but not the joint
 * allele fraction / copy ratio segmenter
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public abstract class ScalarHMMSegmenter<DATA> extends ClusteringGenomicHMMSegmenter<DATA, Double> {

    private static final double DEFAULT_INITIAL_CONCENTRATION = 1;

    protected static final int NEUTRAL_VALUE_INDEX = 0;
    private static final double DEFAULT_MEMORY_LENGTH = 5e7;

    public ScalarHMMSegmenter(final int initialNumStates, final List<SimpleInterval> positions, final List<DATA> data, final List<Double> initialHiddenStates) {
        super(positions, data, initialHiddenStates, uniformWeights(initialNumStates), DEFAULT_INITIAL_CONCENTRATION, DEFAULT_MEMORY_LENGTH);
    }

    @Override
    protected boolean hiddenStateValuesHaveConverged(final List<Double> oldHiddenStateValues) {
        final double[] oldHiddenStateArray = Doubles.toArray(oldHiddenStateValues);
        final double[] newHiddenStateArray = Doubles.toArray(hiddenStateValues);
        return oldHiddenStateArray.length == newHiddenStateArray.length &&
                GATKProtectedMathUtils.maxDifference(oldHiddenStateArray, newHiddenStateArray) < CONVERGENCE_THRESHOLD;
    }

    // filter out components that have low weight and are too close to another component -- these will
    // die out eventually in EM, but very slowly, so we hasten their demise for quicker convergence
    @Override
    protected void pruneUnusedComponents() {
        final int K = weights.length;

        final Set<Integer> componentsToPrune = new TreeSet<>();
        for (final int state : IntStream.range(0, K).filter(n -> n != NEUTRAL_VALUE_INDEX).toArray()) {
            final int closestOtherState = MathUtils.maxElementIndex(IntStream.range(0, K)
                    .mapToDouble(n -> n == state ? Double.NEGATIVE_INFINITY : -Math.abs(hiddenStateValues.get(n) - hiddenStateValues.get(state)))
                    .toArray());
            final boolean hasLowWeight = weights[state] < MAX_WEIGHT_CONSIDERED_FOR_PRUNING;
            final boolean isCloseToNeighbor = Math.abs(hiddenStateValues.get(state) - hiddenStateValues.get(closestOtherState)) < DISTANCE_TO_NEIGHBOR_TO_BE_CONSIDERED_SPURIOUS;
            final boolean isVeryCloseToNeighbor = Math.abs(hiddenStateValues.get(state) - hiddenStateValues.get(closestOtherState)) < DISTANCE_TO_NEIGHBOR_TO_BE_CONSIDERED_DEFINITELY_SPURIOUS;

            final boolean hasLessWeightThanNeighbor = weights[state] < weights[closestOtherState];

            final boolean tinyWeight = weights[state] < AUTOMATICALLY_PRUNED_WEIGHT;
            final boolean parasite = hasLowWeight && isCloseToNeighbor && hasLessWeightThanNeighbor;
            final boolean isClone = isVeryCloseToNeighbor && hasLessWeightThanNeighbor;
            if ( tinyWeight || parasite || isClone) {
                componentsToPrune.add(state);
            }
        }

        weights = IntStream.range(0, K)
                .filter(n -> !componentsToPrune.contains(n)).mapToDouble(n -> weights[n]).toArray();
        hiddenStateValues = IntStream.range(0, K)
                .filter(n -> !componentsToPrune.contains(n)).mapToObj(hiddenStateValues::get).collect(Collectors.toList());

        final StringBuilder message = new StringBuilder("After pruning (state, weight) pairs are: ");
        for (int n = 0; n < weights.length; n++) {
            final double weight = weights[n];
            message.append(String.format("(%f, %f)", hiddenStateValues.get(n), weights[n]) + ((n < weights.length - 1) ? "; " : "."));
        }
        logger.info(message);
    }

    private static double[] uniformWeights(final int numStates) {
        return Collections.nCopies(numStates, 1.0/numStates).stream().mapToDouble(x->x).toArray();
    }

    @Override
    protected void relearnHiddenStateValues(final ExpectationStep eStep) {
        final ClusteringGenomicHMM<DATA, Double> model = makeModel();
        // by convention, state = 0 represents the neutral value (minor allele fraction = 1/2 or copy ratio = 1)
        // which we always retain and do not wish to adjust via MLE.  Thus we start at state 1
        for (final int state : IntStream.range(0, hiddenStateValues.size()).filter(n -> n != NEUTRAL_VALUE_INDEX).toArray()) {
            final Function<Double, Double> objective = f -> IntStream.range(0, data.size())
                    .filter(n -> eStep.pStateAtPosition(state, n) > NEGLIGIBLE_POSTERIOR_FOR_M_STEP)
                    .mapToDouble(n -> eStep.pStateAtPosition(state, n) * model.logEmissionProbability(data.get(n), f))
                    .sum();
            hiddenStateValues.set(state, OptimizationUtils.singleNewtonArgmaxUpdate(objective, minHiddenStateValue(),
                    maxHiddenStateValue(), hiddenStateValues.get(state)));
        }
    }

    protected abstract double minHiddenStateValue();
    protected abstract double maxHiddenStateValue();
}
