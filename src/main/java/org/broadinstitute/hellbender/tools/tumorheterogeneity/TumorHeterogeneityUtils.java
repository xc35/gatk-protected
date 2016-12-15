package org.broadinstitute.hellbender.tools.tumorheterogeneity;

import com.google.common.collect.Sets;
import org.apache.commons.collections4.ListUtils;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.special.Gamma;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.PopulationMixture.VariantProfileCollection;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyState;
import org.broadinstitute.hellbender.tools.tumorheterogeneity.ploidystate.PloidyStatePrior;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.mcmc.coordinates.CoordinateUtils;
import org.broadinstitute.hellbender.utils.mcmc.coordinates.SimplexPosition;
import org.broadinstitute.hellbender.utils.mcmc.coordinates.WalkerPosition;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Helper class containing package-private methods used in org.broadinstitute.hellbender.tools.tumorheterogeneity.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class TumorHeterogeneityUtils {
    private static final Logger logger = LogManager.getLogger(TumorHeterogeneityUtils.class);

    static final int NUM_POPULATIONS_CLONAL = 2;

    static final double EPSILON = 1E-10;
    static final double COPY_RATIO_EPSILON = 1E-3; //below this, use mirrored minor-allele fraction posterior

    static final double CONCENTRATION_MIN = EPSILON;
    static final double CONCENTRATION_MAX = 1. + EPSILON;

    static final double COPY_RATIO_NOISE_CONSTANT_MIN = EPSILON;
    static final double COPY_RATIO_NOISE_CONSTANT_MAX = 5E-2;

    static final double COPY_RATIO_NOISE_FACTOR_MIN = EPSILON;
    static final double COPY_RATIO_NOISE_FACTOR_MAX = 1.;

    static final double MINOR_ALLELE_FRACTION_NOISE_FACTOR_MIN = EPSILON;
    static final double MINOR_ALLELE_FRACTION_NOISE_FACTOR_MAX = 1.;

    static final double PLOIDY_MIN = EPSILON;

    private static final int CONCENTRATION_WALKER_DIMENSION_INDEX = 0;
    private static final int COPY_RATIO_NOISE_CONSTANT_WALKER_DIMENSION_INDEX = 1;
    private static final int COPY_RATIO_NOISE_FACTOR_WALKER_DIMENSION_INDEX = 2;
    private static final int MINOR_ALLELE_FRACTION_NOISE_FACTOR_WALKER_DIMENSION_INDEX = 3;
    private static final int INITIAL_PLOIDY_WALKER_DIMENSION_INDEX = 4;
    private static final int POPULATION_FRACTIONS_WALKER_DIMENSION_START_INDEX = 5;
    static final int NUM_GLOBAL_PARAMETERS = 5;

    private TumorHeterogeneityUtils() {}

    static double calculateLogPosterior(final TumorHeterogeneityState state,
                                        final TumorHeterogeneityData data) {
        if (isOutsideBounds(state, data)) {
            return Double.NEGATIVE_INFINITY;
        }

        //copy-ratio noise-constant prior
        final double copyRatioNoiseConstantPriorAlpha = data.priors().copyRatioNoiseConstantPriorHyperparameterValues().getAlpha();
        final double copyRatioNoiseConstantPriorBeta = data.priors().copyRatioNoiseConstantPriorHyperparameterValues().getBeta();
        final double copyRatioNoiseConstant = state.copyRatioNoiseConstant();
        final double logPriorCopyRatioNoiseConstant =
                copyRatioNoiseConstantPriorAlpha * Math.log(copyRatioNoiseConstantPriorBeta)
                        + (copyRatioNoiseConstantPriorAlpha - 1.) * Math.log(copyRatioNoiseConstant)
                        - copyRatioNoiseConstantPriorBeta * copyRatioNoiseConstant
                        - Gamma.logGamma(copyRatioNoiseConstantPriorAlpha);

        //copy-ratio noise-factor prior
        final double copyRatioNoiseFactorPriorAlpha = data.priors().copyRatioNoiseFactorPriorHyperparameterValues().getAlpha();
        final double copyRatioNoiseFactorPriorBeta = data.priors().copyRatioNoiseFactorPriorHyperparameterValues().getBeta();
        final double copyRatioNoiseFactor = state.copyRatioNoiseFactor();
        final double logPriorCopyRatioNoiseFactor =
                (copyRatioNoiseFactorPriorAlpha - 1.) * Math.log(Math.max(EPSILON, copyRatioNoiseFactor))
                        + (copyRatioNoiseFactorPriorBeta - 1.) * Math.log(Math.max(EPSILON, 1. - copyRatioNoiseFactor));

        //minor-allele-fraction noise-factor prior
        final double minorAlleleFractionNoiseFactorPriorAlpha = data.priors().minorAlleleFractionNoiseFactorPriorHyperparameterValues().getAlpha();
        final double minorAlleleFractionNoiseFactorPriorBeta = data.priors().minorAlleleFractionNoiseFactorPriorHyperparameterValues().getBeta();
        final double minorAlleleFractionNoiseFactor = state.minorAlleleFractionNoiseFactor();
        final double logPriorMinorAlleleFractionNoiseFactor =
                (minorAlleleFractionNoiseFactorPriorAlpha - 1.) * Math.log(Math.max(EPSILON, minorAlleleFractionNoiseFactor))
                        + (minorAlleleFractionNoiseFactorPriorBeta - 1.) * Math.log(Math.max(EPSILON, 1. - minorAlleleFractionNoiseFactor));

        //variant-profiles prior
        final VariantProfileCollection variantProfileCollection = state.populationMixture().variantProfileCollection();
        double logPriorVariantProfiles = 0.;
        for (int populationIndex = 0; populationIndex < variantProfileCollection.numVariantPopulations(); populationIndex++) {
            for (int segmentIndex = 0; segmentIndex < variantProfileCollection.numSegments(); segmentIndex++) {
                final PloidyState ploidyState = variantProfileCollection.ploidyState(populationIndex, segmentIndex);
                logPriorVariantProfiles += data.length(segmentIndex) * data.priors().ploidyStatePrior().logProbability(ploidyState);
            }
        }

        //copy-ratio--minor-allele-fraction likelihood
        double logLikelihoodSegments = 0.;
        final double ploidy = state.ploidy();
        for (int segmentIndex = 0; segmentIndex < data.numSegments(); segmentIndex++) {
            final double mAlleleCopyNumber = state.populationMixture().calculatePopulationAveragedCopyNumberFunction(segmentIndex, PloidyState::m);
            final double nAlleleCopyNumber = state.populationMixture().calculatePopulationAveragedCopyNumberFunction(segmentIndex, PloidyState::n);
            final double totalCopyNumber = mAlleleCopyNumber + nAlleleCopyNumber;
            final double copyRatio = totalCopyNumber / Math.max(EPSILON, ploidy);
            final double minorAlleleFraction = calculateMinorAlleleFraction(mAlleleCopyNumber, nAlleleCopyNumber);
            logLikelihoodSegments += data.logDensity(
                    segmentIndex, copyRatio, minorAlleleFraction,
                    state.copyRatioNoiseConstant(), state.copyRatioNoiseFactor(), state.minorAlleleFractionNoiseFactor());
        }

        //log penalty for mismatch between initial ploidy and ploidy resulting from proposeVariantProfileCollection
        double logPloidyMismatchPenalty = -data.priors().ploidyMismatchPenalty() * Math.abs(state.initialPloidy() - state.ploidy());

        logger.debug("Log-posterior components:"
                + " " + logPriorCopyRatioNoiseConstant
                + " " + logPriorCopyRatioNoiseFactor
                + " " + logPriorMinorAlleleFractionNoiseFactor
                + " " + logPriorVariantProfiles
                + " " + logLikelihoodSegments
                + " " + logPloidyMismatchPenalty);

        return logPriorCopyRatioNoiseConstant + logPriorCopyRatioNoiseFactor  + logPriorMinorAlleleFractionNoiseFactor
                + logPriorVariantProfiles + logLikelihoodSegments + logPloidyMismatchPenalty;
    }

    /**
     * Calculates the log of the Jacobian factor for the state-to-walker transformation.
     */
    static double calculateLogJacobianFactor(final TumorHeterogeneityState state,
                                             final TumorHeterogeneityData data) {
        final double ploidyMax = data.priors().ploidyStatePrior().maxCopyNumber();
        return CoordinateUtils.calculateLogJacobianFactor(state.concentration(), CONCENTRATION_MIN, CONCENTRATION_MAX)
                + CoordinateUtils.calculateLogJacobianFactor(state.copyRatioNoiseConstant(), COPY_RATIO_NOISE_CONSTANT_MIN, COPY_RATIO_NOISE_CONSTANT_MAX)
                + CoordinateUtils.calculateLogJacobianFactor(state.copyRatioNoiseFactor(), COPY_RATIO_NOISE_FACTOR_MIN, COPY_RATIO_NOISE_FACTOR_MAX)
                + CoordinateUtils.calculateLogJacobianFactor(state.minorAlleleFractionNoiseFactor(), MINOR_ALLELE_FRACTION_NOISE_FACTOR_MIN, MINOR_ALLELE_FRACTION_NOISE_FACTOR_MAX)
                + CoordinateUtils.calculateLogJacobianFactor(state.populationMixture().ploidy(data), PLOIDY_MIN, ploidyMax)
                + SimplexPosition.calculateLogJacobianFactor(state.populationMixture().populationFractions());
    }

    /**
     * Transforms a walker position (i.e., a point in unbounded N-dimensional space) to
     * a {@link TumorHeterogeneityState} (which is composed of parameters that may be bounded or confined to the
     * unit simplex).  This transformation includes a proposal step to randomly sample a {@link VariantProfileCollection}
     * conditional on the other model parameters and a proposed "initial" ploidy
     * using {@link TumorHeterogeneityUtils#proposeVariantProfileCollection}, so the transformation is not
     * deterministic.
     */
    static TumorHeterogeneityState transformWalkerPositionToState(final WalkerPosition walkerPosition,
                                                                  final RandomGenerator rng,
                                                                  final TumorHeterogeneityData data,
                                                                  final List<List<Integer>> totalCopyNumberProductStates,
                                                                  final Map<Integer, Set<PloidyState>> ploidyStateSetsMap) {
        final double concentration = CoordinateUtils.transformWalkerCoordinateToBoundedVariable(
                walkerPosition.get(CONCENTRATION_WALKER_DIMENSION_INDEX), CONCENTRATION_MIN, CONCENTRATION_MAX);
        final double copyRatioNoiseConstant = CoordinateUtils.transformWalkerCoordinateToBoundedVariable(
                walkerPosition.get(COPY_RATIO_NOISE_CONSTANT_WALKER_DIMENSION_INDEX), COPY_RATIO_NOISE_CONSTANT_MIN, COPY_RATIO_NOISE_CONSTANT_MAX);
        final double copyRatioNoiseFactor = CoordinateUtils.transformWalkerCoordinateToBoundedVariable(
                walkerPosition.get(COPY_RATIO_NOISE_FACTOR_WALKER_DIMENSION_INDEX), COPY_RATIO_NOISE_FACTOR_MIN, COPY_RATIO_NOISE_FACTOR_MAX);
        final double minorAlleleFractionNoiseFactor = CoordinateUtils.transformWalkerCoordinateToBoundedVariable(
                walkerPosition.get(MINOR_ALLELE_FRACTION_NOISE_FACTOR_WALKER_DIMENSION_INDEX), MINOR_ALLELE_FRACTION_NOISE_FACTOR_MIN, MINOR_ALLELE_FRACTION_NOISE_FACTOR_MAX);
        final double ploidyMax = data.priors().ploidyStatePrior().maxCopyNumber();
        final double initialPloidy = CoordinateUtils.transformWalkerCoordinateToBoundedVariable(
                walkerPosition.get(INITIAL_PLOIDY_WALKER_DIMENSION_INDEX), PLOIDY_MIN, ploidyMax);

        final PopulationMixture.PopulationFractions populationFractions = new PopulationMixture.PopulationFractions(
                SimplexPosition.calculateSimplexPositionFromWalkerPosition(
                        new WalkerPosition(walkerPosition.subList(POPULATION_FRACTIONS_WALKER_DIMENSION_START_INDEX, walkerPosition.numDimensions()))));

        final VariantProfileCollection variantProfileCollection = proposeVariantProfileCollection(
                rng, copyRatioNoiseConstant, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, initialPloidy,
                populationFractions, data, totalCopyNumberProductStates, ploidyStateSetsMap);

        final PopulationMixture populationMixture = new PopulationMixture(populationFractions, variantProfileCollection, data.priors().normalPloidyState());
        final double ploidy = populationMixture.ploidy(data);

        return new TumorHeterogeneityState(concentration, copyRatioNoiseConstant, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor, initialPloidy, ploidy, populationMixture);
    }

    /**
     * Transforms a {@link TumorHeterogeneityState} to a walker position.  This is only used to initialize a ball
     * of walkers around an initial state in {@link TumorHeterogeneityModeller#initializeWalkerBall}.
     */
    static WalkerPosition transformStateToWalkerPosition(final TumorHeterogeneityState state,
                                                         final TumorHeterogeneityData data) {
        final double concentrationWalkerCoordinate = CoordinateUtils.transformBoundedVariableToWalkerCoordinate(
                state.concentration(), CONCENTRATION_MIN, CONCENTRATION_MAX);
        final double copyRatioNoiseConstantWalkerCoordinate = CoordinateUtils.transformBoundedVariableToWalkerCoordinate(
                state.copyRatioNoiseConstant(), COPY_RATIO_NOISE_CONSTANT_MIN, COPY_RATIO_NOISE_CONSTANT_MAX);
        final double copyRatioNoiseFactorWalkerCoordinate = CoordinateUtils.transformBoundedVariableToWalkerCoordinate(
                state.copyRatioNoiseFactor(), COPY_RATIO_NOISE_FACTOR_MIN, COPY_RATIO_NOISE_FACTOR_MAX);
        final double minorAlleleFractionNoiseFactorWalkerCoordinate = CoordinateUtils.transformBoundedVariableToWalkerCoordinate(
                state.minorAlleleFractionNoiseFactor(), MINOR_ALLELE_FRACTION_NOISE_FACTOR_MIN, MINOR_ALLELE_FRACTION_NOISE_FACTOR_MAX);
        final double ploidyMax = data.priors().ploidyStatePrior().maxCopyNumber();
        final double initialPloidyWalkerCoordinate = CoordinateUtils.transformBoundedVariableToWalkerCoordinate(
                state.initialPloidy(), PLOIDY_MIN, ploidyMax);
        final WalkerPosition populationFractionsWalkerCoordinates = SimplexPosition.calculateWalkerPositionFromSimplexPosition(state.populationMixture().populationFractions());

        return new WalkerPosition(ListUtils.union(
                Arrays.asList(
                        concentrationWalkerCoordinate,
                        copyRatioNoiseConstantWalkerCoordinate,
                        copyRatioNoiseFactorWalkerCoordinate,
                        minorAlleleFractionNoiseFactorWalkerCoordinate,
                        initialPloidyWalkerCoordinate),
                populationFractionsWalkerCoordinates));
    }

    /**
     * Conditional on the global model parameters, population fractions, and a proposed "initial" ploidy, heuristically
     * samples a {@link VariantProfileCollection} from the posterior distribution.
     */
    private static VariantProfileCollection proposeVariantProfileCollection(final RandomGenerator rng,
                                                                            final double copyRatioNoiseConstant,
                                                                            final double copyRatioNoiseFactor,
                                                                            final double minorAlleleFractionNoiseFactor,
                                                                            final double initialPloidy,
                                                                            final PopulationMixture.PopulationFractions populationFractions,
                                                                            final TumorHeterogeneityData data,
                                                                            final List<List<Integer>> totalCopyNumberProductStates,
                                                                            final Map<Integer, Set<PloidyState>> ploidyStateSetsMap) {
        final int numVariantPopulations = populationFractions.numVariantPopulations();
        final int numSegments = data.numSegments();
        final PloidyState normalPloidyState = data.priors().normalPloidyState();
        final PloidyStatePrior ploidyStatePrior = data.priors().ploidyStatePrior();

        //initialize a list of variant profiles to store the result
        final List<PopulationMixture.VariantProfile> variantProfiles =
                TumorHeterogeneityState.initializeNormalProfiles(numVariantPopulations, numSegments, normalPloidyState);

        for (int segmentIndex = 0; segmentIndex < numSegments; segmentIndex++) {
            final int si = segmentIndex;

            //for all possible copy-number product states, calculate the copy-ratio likelihood given the proposed ploidy
            final List<Double> logProbabilitiesCopyRatio = totalCopyNumberProductStates.stream()
                    .map(tcnps -> calculateTotalCopyNumber(populationFractions, tcnps, normalPloidyState) / Math.max(EPSILON, initialPloidy))
                    .map(cr -> data.copyRatioLogDensity(si, cr, copyRatioNoiseConstant, copyRatioNoiseFactor))
                    .collect(Collectors.toList());
            //take maximum likelihood copy-number product states according to their copy-ratio--only likelihoods
            final int maxLikelihoodCopyNumberProductStateIndex = IntStream.range(0, totalCopyNumberProductStates.size()).boxed()
                    .max((i, j) -> Double.compare(logProbabilitiesCopyRatio.get(i), logProbabilitiesCopyRatio.get(j))).get();
            final List<Integer> totalCopyNumberProductState = totalCopyNumberProductStates.get(maxLikelihoodCopyNumberProductStateIndex);
            //calculate the copy ratio of the sampled copy-number product state
            final double totalCopyRatio = calculateTotalCopyNumber(populationFractions, totalCopyNumberProductState, normalPloidyState) / Math.max(EPSILON, initialPloidy);

            //for all ploidy-state product states consistent with the sampled copy-number product state, calculate the copy-ratio--minor-allele-fraction posteriors
            final List<List<PloidyState>> ploidyStateProductStates =
                    new ArrayList<>(Sets.cartesianProduct(totalCopyNumberProductState.stream().map(ploidyStateSetsMap::get).collect(Collectors.toList())));
            final List<Double> logProbabilitiesPloidyStateProductStates = ploidyStateProductStates.stream()
                    .map(psps -> data.logDensity(si, totalCopyRatio, calculateMinorAlleleFraction(populationFractions, psps, normalPloidyState), copyRatioNoiseConstant, copyRatioNoiseFactor, minorAlleleFractionNoiseFactor)
                            + psps.stream().mapToDouble(ps -> data.length(si) * ploidyStatePrior.logProbability(ps)).sum())
                    .collect(Collectors.toList());
            //take maximum a posteriori ploidy-state product state according to their copy-ratio--minor-allele-fraction posteriors
            final int maxPosteriorPloidyStateProductStateIndex = IntStream.range(0, ploidyStateProductStates.size()).boxed()
                    .max((i, j) -> Double.compare(logProbabilitiesPloidyStateProductStates.get(i), logProbabilitiesPloidyStateProductStates.get(j))).get();
            final List<PloidyState> ploidyStateProductState = ploidyStateProductStates.get(maxPosteriorPloidyStateProductStateIndex);

            //store the sampled ploidy-state product state in the list of variant profiles
            IntStream.range(0, numVariantPopulations).forEach(i -> variantProfiles.get(i).set(si, ploidyStateProductState.get(i)));
        }
        //return a new VariantProfileCollection created from the list of variant profiles
        return new VariantProfileCollection(variantProfiles);
    }

    private static double calculateTotalCopyNumber(final PopulationMixture.PopulationFractions populationFractions,
                                                   final List<Integer> totalCopyNumberProductState,
                                                   final PloidyState normalPloidyState) {
        final int numVariantPopulations = populationFractions.numVariantPopulations();
        //loop performs better than stream
        double totalCopyNumber = 0.;
        for (int populationIndex = 0; populationIndex < numVariantPopulations; populationIndex++) {
            totalCopyNumber += totalCopyNumberProductState.get(populationIndex) * populationFractions.get(populationIndex);
        }
        totalCopyNumber += normalPloidyState.total() * populationFractions.normalFraction();
        return totalCopyNumber;
    }

    private static double calculateMinorAlleleFraction(final PopulationMixture.PopulationFractions populationFractions,
                                                       final List<PloidyState> ploidyStateProductState,
                                                       final PloidyState normalPloidyState) {
        final int numVariantPopulations = populationFractions.numVariantPopulations();
        //performance is not as critical here, so we use streams
        final double mAlleleCopyNumber = IntStream.range(0, numVariantPopulations)
                .mapToDouble(i -> ploidyStateProductState.get(i).m() * populationFractions.get(i))
                .sum() + normalPloidyState.m() * populationFractions.normalFraction();
        final double nAlleleCopyNumber = IntStream.range(0, numVariantPopulations)
                .mapToDouble(i -> ploidyStateProductState.get(i).n() * populationFractions.get(i))
                .sum() + normalPloidyState.n() * populationFractions.normalFraction();
        return calculateMinorAlleleFraction(mAlleleCopyNumber, nAlleleCopyNumber);
    }

    private static double calculateMinorAlleleFraction(final double m, final double n) {
        return Math.min(m, n) / Math.max(EPSILON, m + n);
    }

    private static boolean isOutsideBounds(final TumorHeterogeneityState state,
                                           final TumorHeterogeneityData data) {
        final double ploidyMax = data.priors().ploidyStatePrior().maxCopyNumber();
        return state.concentration() < CONCENTRATION_MIN || state.concentration() > CONCENTRATION_MAX ||
                state.copyRatioNoiseConstant() < COPY_RATIO_NOISE_CONSTANT_MIN || state.copyRatioNoiseConstant() > COPY_RATIO_NOISE_CONSTANT_MAX ||
                state.copyRatioNoiseFactor() < COPY_RATIO_NOISE_FACTOR_MIN || state.copyRatioNoiseFactor() > COPY_RATIO_NOISE_FACTOR_MAX ||
                state.minorAlleleFractionNoiseFactor() < MINOR_ALLELE_FRACTION_NOISE_FACTOR_MIN || state.minorAlleleFractionNoiseFactor() > MINOR_ALLELE_FRACTION_NOISE_FACTOR_MAX ||
                state.ploidy() < PLOIDY_MIN || state.ploidy() > ploidyMax ||
                state.populationMixture().variantProfileCollection().stream().anyMatch(PopulationMixture.VariantProfile::isCompleteDeletion);
    }
}