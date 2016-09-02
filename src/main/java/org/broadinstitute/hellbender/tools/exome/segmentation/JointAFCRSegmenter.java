package org.broadinstitute.hellbender.tools.exome.segmentation;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.primitives.Doubles;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.ml.clustering.DBSCANClusterer;
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionGlobalParameters;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by davidben on 9/6/16.
 */
public class JointAFCRSegmenter extends ClusteringGenomicHMMSegmenter<TargetOrHetDatum, AFCRHiddenState> {
    final AlleleFractionGlobalParameters parameters;
    final AllelicPanelOfNormals allelicPoN;
    final double log2CoverageCauchyWidth;

    // make this a parameter
    private static final int MAX_NUM_STATES = 50;

    private static final Comparator<Locatable> INTERVAL_COMPARATOR = IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR;

    private JointAFCRSegmenter(final List<SimpleInterval> positions,
                              final List<TargetOrHetDatum> data,
                              final List<AFCRHiddenState> hiddenStateValues,
                              final double[] weights,
                              final double concentration,
                              final double memoryLength,
                              final AlleleFractionGlobalParameters parameters,
                              final AllelicPanelOfNormals allelicPoN,
                              final double log2CoverageCauchyWidth) {

        super(positions, data, hiddenStateValues, weights, concentration, memoryLength);
        this.parameters = parameters;
        this.allelicPoN = allelicPoN;
        this.log2CoverageCauchyWidth = log2CoverageCauchyWidth;
        logStatesAndWeights();
    }

    private void logStatesAndWeights() {
        final StringBuilder message = new StringBuilder("Joint segmenter has following (copy-ratio, allele-fraction): weight pairs: ");
        for (int n = 0; n < weights.length; n++) {
            final AFCRHiddenState state = hiddenStateValues.get(n);
            final double cr = state.getLog2CopyRatio();
            final double af = state.getMinorAlleleFraction();
            final double weight = weights[n];
            message.append(String.format("(%f, %f): %f", cr, af, weight) + ((n < weights.length - 1) ? ";" : "."));
        }
        logger.info(message);
    }

    public static JointAFCRSegmenter createJointSegmenter(final CopyRatioSegmenter copyRatioSegmenter, final AlleleFractionSegmenter alleleFractionSegmenter) {
        copyRatioSegmenter.makeSureParametersHaveBeenLearned();
        alleleFractionSegmenter.makeSureParametersHaveBeenLearned();
        List<IntervalDatumIndex> targetsWithCopyRatio = IntStream.range(0, copyRatioSegmenter.numPositions())
                .mapToObj(n -> new IntervalDatumIndex(copyRatioSegmenter.getPosition(n), new TargetOrHetDatum(copyRatioSegmenter.getDatum(n)), n))
                .collect(Collectors.toList());
        List<IntervalDatumIndex> hetsWithAllelicCount = IntStream.range(0, alleleFractionSegmenter.numPositions())
                .mapToObj(n -> new IntervalDatumIndex(alleleFractionSegmenter.getPosition(n), new TargetOrHetDatum(alleleFractionSegmenter.getDatum(n)), n))
                .collect(Collectors.toList());
        List<IntervalDatumIndex> positionsWithData = CollectionUtils.collate(targetsWithCopyRatio, hetsWithAllelicCount);

        final List<SimpleInterval> allPositions = positionsWithData.stream().map(IntervalDatumIndex::getInterval).collect(Collectors.toList());
        final List<TargetOrHetDatum> allData = positionsWithData.stream().map(IntervalDatumIndex::getDatum).collect(Collectors.toList());


        // pairs of targets with immediately preceding and following hets, and vice versa
        final List<ImmutablePair<Integer, Integer>> neighboringTargetHetIndices = new ArrayList<>(2 * allData.size());
        final int numHets = alleleFractionSegmenter.numPositions();
        final int numTargets = copyRatioSegmenter.numPositions();

        for (int n = 0; n < allData.size(); n++) {
            if (allData.get(n).isTarget()) {
                final int targetIndex = positionsWithData.get(n).index;
                final int followingHetIndex = n - targetIndex;
                final int precedingHetIndex = followingHetIndex - 1;
                neighboringTargetHetIndices.add(new ImmutablePair<>(targetIndex, precedingHetIndex));
                neighboringTargetHetIndices.add(new ImmutablePair<>(targetIndex, followingHetIndex));
            } else {
                final int hetIndex = positionsWithData.get(n).index;
                final int followingTargetIndex = n - hetIndex;
                final int precedingTargetIndex = followingTargetIndex - 1;
                neighboringTargetHetIndices.add(new ImmutablePair<>(precedingTargetIndex, hetIndex));
                neighboringTargetHetIndices.add(new ImmutablePair<>(followingTargetIndex, hetIndex));
            }
        }

        // fill a matrix (# of CR hidden state x # of AF hidden states) of co-occuring CR and AF states
        final CopyRatioSegmenter.ExpectationStep copyRatioEStep = copyRatioSegmenter.getExpectationStep();
        final AlleleFractionSegmenter.ExpectationStep alleleFractionEStep = alleleFractionSegmenter.getExpectationStep();

        final int numCopyRatioStates = copyRatioSegmenter.hiddenStateValues.size();
        final int numAlleleFractionStates = alleleFractionSegmenter.hiddenStateValues.size();
        final double[][] cooccurrence = new double[copyRatioSegmenter.hiddenStateValues.size()][alleleFractionSegmenter.hiddenStateValues.size()];

        neighboringTargetHetIndices.stream()
                .filter(p -> p.getLeft() >= 0 && p.getRight() >= 0 && p.getLeft() < numTargets && p.getRight() < numHets)
                .filter(p -> copyRatioSegmenter.positions.get(p.getLeft()).getContig().equals(alleleFractionSegmenter.positions.get(p.getRight()).getContig()))
                .forEach(p -> {
                    final int target = p.getLeft();
                    final int het = p.getRight();
                    for (int i = 0; i < numCopyRatioStates; i++) {
                        for (int j = 0; j < numAlleleFractionStates; j++) {
                            cooccurrence[i][j] += copyRatioEStep.pStateAtPosition(i, target) * alleleFractionEStep.pStateAtPosition(j, het);
                        }
                    }
                });

        final List<ImmutablePair<AFCRHiddenState, Double>> jointStatesAndWeights = new ArrayList<>();
        for (int i = 0; i < numCopyRatioStates; i++) {
            for (int j = 0; j < numAlleleFractionStates; j++) {
                jointStatesAndWeights.add(new ImmutablePair<>(
                        new AFCRHiddenState(alleleFractionSegmenter.hiddenStateValues.get(j), copyRatioSegmenter.hiddenStateValues.get(i)), cooccurrence[i][j]));
            }
        }

        // descending order by cooccurrence score
        Collections.sort(jointStatesAndWeights, (jsaw1, jsaw2) -> Doubles.compare(jsaw2.getRight(), jsaw1.getRight()));
        final List<ImmutablePair<AFCRHiddenState, Double>> bestJointStatesAndWeights =
                jointStatesAndWeights.subList(0, Math.min(MAX_NUM_STATES, jointStatesAndWeights.size()));
        final List<AFCRHiddenState> hiddenStates = bestJointStatesAndWeights.stream()
                .map(jsaw -> jsaw.getLeft()).collect(Collectors.toList());
        final double[] unnormalizedWeights = bestJointStatesAndWeights.stream()
                .mapToDouble(jsaw -> jsaw.getRight()).toArray();
        final double[] weights = MathUtils.normalizeFromRealSpace(unnormalizedWeights);
        final double concentration = geometricMean(copyRatioSegmenter.getConcentration(), alleleFractionSegmenter.getConcentration());
        final double memoryLength = geometricMean(copyRatioSegmenter.getMemoryLength(), alleleFractionSegmenter.getMemoryLength());
        return new JointAFCRSegmenter(allPositions, allData, hiddenStates, weights, concentration, memoryLength,
                alleleFractionSegmenter.getBiasParameters(), alleleFractionSegmenter.getAllelicPoN(), copyRatioSegmenter.getLogCoverageCauchyWidth());
    }

    // joint segmenter does not relearn any parameters
    @Override
    protected void relearnAdditionalParameters(final ExpectationStep eStep) { }

    @Override
    protected void relearnHiddenStateValues(final ExpectationStep eStep) { }

    @Override
    protected boolean hiddenStateValuesHaveConverged(final List<AFCRHiddenState> hiddenStateValues) { return true; }


    //TODO: consider non-trivial pruning based on weights
    @Override
    protected void pruneUnusedComponents() {
        logStatesAndWeights();
        //DBSCANClusterer clusterer = new DBSCANClusterer()
    }

    @Override
    protected ClusteringGenomicHMM<TargetOrHetDatum, AFCRHiddenState> makeModel() {
        return new JointAFCRHiddenMarkovModel(hiddenStateValues, weights, getMemoryLength(), parameters, allelicPoN, log2CoverageCauchyWidth);
    }

    // store an interval, a target/het datum, and the index of that datum within either the targets or hets,
    // depending on the type of datum
    private static final class IntervalDatumIndex implements Comparable<IntervalDatumIndex> {
        private final SimpleInterval interval;
        private final TargetOrHetDatum datum;
        public final int index;

        public IntervalDatumIndex(SimpleInterval interval, TargetOrHetDatum datum, int index) {
            this.interval = interval;
            this.datum = datum;
            this.index = index;
        }

        @Override
        public int compareTo(final IntervalDatumIndex other) {
            return INTERVAL_COMPARATOR.compare(this.interval, other.interval);
        }

        public SimpleInterval getInterval() { return interval; }
        public TargetOrHetDatum getDatum() { return  datum; }
    }

    private static double geometricMean(final double x, final double y) {
        return Math.sqrt(x) * Math.sqrt(y);
    }

}
