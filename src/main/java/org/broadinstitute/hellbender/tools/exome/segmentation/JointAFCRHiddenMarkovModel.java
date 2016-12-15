package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionGlobalParameters;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public class JointAFCRHiddenMarkovModel extends ClusteringGenomicHMM<TargetOrHetDatum, AFCRHiddenState> {
    private final AlleleFractionGlobalParameters parameters;
    private final AllelicPanelOfNormals allelicPoN;
    private final double logCoverageCauchyWidth;

    public JointAFCRHiddenMarkovModel(final List<AFCRHiddenState> hiddenStateValues, final double[] weights,
                                      final double memoryLength, final AlleleFractionGlobalParameters parameters,
                                      final AllelicPanelOfNormals allelicPoN, final double logCoverageCauchyWidth) {
        super(hiddenStateValues, weights, memoryLength);
        this.parameters = parameters;
        this.allelicPoN = allelicPoN;
        this.logCoverageCauchyWidth = logCoverageCauchyWidth;
    }

    @Override
    public double logEmissionProbability(final TargetOrHetDatum datum, final Integer state, final SimpleInterval position) {
        return logEmissionProbability(datum, hiddenStateValues.get(state));
    }

    @Override
    public double logEmissionProbability(final TargetOrHetDatum datum, final AFCRHiddenState hiddenState) {
        return datum.isTarget() ?
                CopyRatioHiddenMarkovModel.logEmissionProbability(datum.getCopyRatio(), hiddenState.getLog2CopyRatio(), logCoverageCauchyWidth)
                : AlleleFractionHiddenMarkovModel.logEmissionProbability(datum.getAllelicCount(), hiddenState.getMinorAlleleFraction(), parameters, allelicPoN);
    }
}
