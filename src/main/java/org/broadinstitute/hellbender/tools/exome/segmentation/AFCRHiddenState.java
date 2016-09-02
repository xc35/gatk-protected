package org.broadinstitute.hellbender.tools.exome.segmentation;

/**
 * Created by davidben on 9/6/16.
 */
public class AFCRHiddenState {
    private final double minorAlleleFraction;
    private final double log2CopyRatio;

    public AFCRHiddenState(double minorAlleleFraction, double log2CopyRatio) {
        this.minorAlleleFraction = minorAlleleFraction;
        this.log2CopyRatio = log2CopyRatio;
    }

    public double getMinorAlleleFraction() {
        return minorAlleleFraction;
    }

    public double getLog2CopyRatio() {
        return log2CopyRatio;
    }
}
