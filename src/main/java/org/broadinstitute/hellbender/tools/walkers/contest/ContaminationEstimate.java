package org.broadinstitute.hellbender.tools.walkers.contest;

import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

/**
 * a class that estimates and stores the contamination values for a site.
 */
final class ContaminationEstimate {
    private double[] log10Likelihoods;   // log10Likelihoods at discrete contamination levels

    // precalculate the 128 possible values of epsilon
    private final static double[] linearSpaceErrorProbs = new IndexRange(0, Byte.MAX_VALUE + 1).mapToDouble(n -> Math.pow(10.0, -n/10.0));

    /**
     * create the contamination estimate, given:
     * @param precision the precision value, to what level are we calculating the contamination
     */
    public ContaminationEstimate(final double precision,
                                 final double maf,
                                 final ReadPileup pileup,
                                 final byte uncontaminatedAllele,
                                 final byte hapmapAlt) {
        Utils.validateArg(0.0 <= maf && maf <= 1.0, () -> "Invalid allele Freq: must be between 0 and 1 (inclusive), maf was " + maf);

        //TODO: this is log, not log10 -- everything is inconsistent!!!!!
        log10Likelihoods = new double[(int)Math.ceil(100/precision)+1];

        final byte[] quals = pileup.getBaseQuals();
        PileupElement
        int qualOffset = 0;
        for (byte base : pileup.getBaseQuals()) {
            final double pError = linearSpaceErrorProbs[quals[qualOffset++]];

            for (int index = 0; index < log10Likelihoods.length; index++) {
                final double contaminationRate = (1.0 - index / (double) log10Likelihoods.length);

                if (base == uncontaminatedAllele) {
                    log10Likelihoods[index] += Math.log((1.0 - contaminationRate) * (1.0 - pError) +
                            contaminationRate * (maf * (1.0 - pError) + (1.0 - maf) * pError/3.0));
                } else if(hapmapAlt == base) {
                    log10Likelihoods[index] += Math.log((1.0 - contaminationRate) * pError / 3.0 +
                            contaminationRate * (maf * pError/3.0 + (1.0 - maf) * (1.0 - pError)));
                }
            }
        }
    }

    public void addSiteData(ContaminationEstimate other) {
        if (other == null) return;
        this.log10Likelihoods = MathArrays.ebeAdd(this.log10Likelihoods, other.log10Likelihoods);
    }
}
