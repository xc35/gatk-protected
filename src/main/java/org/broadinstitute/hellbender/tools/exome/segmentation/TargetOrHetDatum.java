package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;

import java.util.Optional;

/**
 * Holds either the {@link org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount} data
 * for a single het site or a double representing a copy ratio
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public class TargetOrHetDatum {
    private final Optional<Double> copyRatio;
    private final Optional<AllelicCount> allelicCount;

    // construct a copy ratio datum
    public TargetOrHetDatum(final double copyRatio) {
        this.copyRatio = Optional.of(copyRatio);
        allelicCount = Optional.empty();
    }

    // construct an allelic count datum
    public TargetOrHetDatum(final AllelicCount allelicCount) {
        this.allelicCount = Optional.of(allelicCount);
        copyRatio = Optional.empty();
    }

    public boolean isTarget() {
        return copyRatio.isPresent();
    }

    public boolean isHet() {
        return allelicCount.isPresent();
    }

    public double getCopyRatio() { return copyRatio.get(); }

    public AllelicCount getAllelicCount() { return allelicCount.get(); }
}
