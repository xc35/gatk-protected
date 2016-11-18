package org.broadinstitute.hellbender.tools.walkers.contest;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

/**
 * a class that estimates and stores the contamination values for a site.
 */
final class ContaminationEstimate {
    private double[] log10Likelihoods;   // log10Likelihoods at discrete contamination levels

    private long basesFor = 0L;
    private long basesAgainst = 0L;
    private final Nucleotide.Counter alleleBreakdown;

    // precalculate the 128 possible values of epsilon
    private final static double[] linearSpaceErrorProbs = new IndexRange(0, Byte.MAX_VALUE + 1).mapToDouble(n -> Math.pow(10.0, -n/10.0));

    /**
     * create the contamination estimate, given:
     * @param precision the precision value, to what level are we calculating the contamination
     */
    public ContaminationEstimate(long basesFor, long basesAgainst, Nucleotide.Counter alleleBreakdown,
                                 final double precision,
                                 final double maf,
                                 final ReadPileup pileup,
                                 final byte uncontaminatedAllele,
                                 final byte hapmapAlt) {

        this.basesFor = basesFor;
        this.basesAgainst = basesAgainst;
        this.alleleBreakdown = alleleBreakdown;

        Utils.validateArg(0.0 <= maf && maf <= 1.0, () -> "Invalid allele Freq: must be between 0 and 1 (inclusive), maf was " + maf);

        //TODO: this is log, not log10 -- everything is inconsistent!!!!!
        log10Likelihoods = new double[(int)Math.ceil(100/precision)+1];

        final byte[] quals = pileup.getBaseQuals();
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

    public double[] getBins() {
        return log10Likelihoods;
    }

    public void add(ContaminationEstimate other) {
        if (other == null) return;
        this.basesFor               += other.basesFor;
        this.basesAgainst           += other.basesAgainst;
        this.alleleBreakdown.increment(other.alleleBreakdown);
        this.log10Likelihoods = MathArrays.ebeAdd(this.log10Likelihoods, other.log10Likelihoods);
    }

    //TODO: there has to be a library to replace this
    public static class ConfidenceInterval {

        private double start;
        private double stop;
        private double contamination;

        public ConfidenceInterval(double log10Probabilities[]) {
            // make a copy of the log10Probabilities in non-log space
            int maxIndex = MathUtils.maxElementIndex(log10Probabilities);

            final double[] probabilities = MathUtils.normalizeFromLog10ToLinearSpace(log10Probabilities);

            double areaUnderCurve = 0;
            int leftIndex = maxIndex;
            int rightIndex = maxIndex;
            while (areaUnderCurve < 0.95) {

                // if the "left" bin is bigger, and can be moved, move it
                if (probabilities[leftIndex] >= probabilities[rightIndex] && leftIndex > 0) {
                    leftIndex--;
                } else {
                    // otherwise move the right bin if possible
                    if (rightIndex < log10Probabilities.length - 1) {
                        rightIndex++;
                    } else {
                        // and if not move the left bin, or die
                        if (leftIndex > 0) {
                            leftIndex--;
                        } else {
                            throw new RuntimeException("Error trying to compute confidence interval");
                        }
                    }
                }

                areaUnderCurve = 0.0;
                for (int x = leftIndex; x <= rightIndex; x++)
                    areaUnderCurve += probabilities[x];
            }
            start = (log10Probabilities.length - rightIndex) * (100.0/log10Probabilities.length);
            stop = (log10Probabilities.length - leftIndex) * (100.0/log10Probabilities.length);
            contamination = (log10Probabilities.length - maxIndex) * (100.0/log10Probabilities.length);
        }

        public double getStart() {
            return start;
        }

        public double getStop() {
            return stop;
        }

        public double getContamination() {
            return contamination;
        }
    }
}
