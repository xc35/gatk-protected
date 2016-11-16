package org.broadinstitute.hellbender.tools.walkers.contest;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * a class that estimates and stores the contamination values for a site.
 */
final class ContaminationEstimate {
    private final double[] bins;                // the bins representing the discrete contamination levels we're evaluating
    private double populationFit = 0.0;
    private String populationName = "";

    // precalculate the 128 possible values of epsilon
    private static double[] precalculatedEpsilon = new IndexRange(0, Byte.MAX_VALUE + 1).mapToDouble(n -> Math.pow(10.0, -n/10.0));

    /**
     * create the contamination estimate, given:
     * @param precision the precision value, to what level are we calculating the contamination
     */
    public ContaminationEstimate(double precision,
                                 final double maf,  //TODO: I think maf == 0 is not an issue but double-check
                                 byte[] bases,
                                 byte[] quals,
                                 byte arrayAllele,
                                 byte hapmapAlt,
                                 String popName,
                                 Locatable locus) {
        Utils.validateArg(0.0 <= maf && maf <= 1.0, () -> "Invalid allele Freq: must be between 0 and 1 (inclusive), maf was " + maf + " for population " + popName);

        bins = new double[(int)Math.ceil(100/precision)+1];
        populationName = popName;
        final double[] realQuals = new IndexRange(0, quals.length).mapToDouble(n -> Math.pow(10.0,-quals[n]/10.0));

        // calculate the contamination for each bin
        int qualOffset = 0;
        for (byte base : bases) {
            final double epsilon = precalculatedEpsilon[quals[qualOffset++]];

            for (int index = 0; index < bins.length; index++) {
                final double contaminationRate = (1.0 - index / (double) bins.length);

                if (base == arrayAllele) {
                    bins[index] += Math.log((1.0 - contaminationRate) * (1.0 - epsilon) +
                            contaminationRate * (maf * (1.0 - epsilon) + (1.0 - maf) * epsilon/3.0));
                    populationFit += Math.log(epsilon);

                } else if(hapmapAlt == base) {
                    bins[index] += Math.log((1.0 - contaminationRate) * epsilon / 3.0 +
                            contaminationRate * (maf * epsilon/3.0 + (1.0 - maf) * (1.0 - epsilon)));

                    populationFit += Math.log(maf + epsilon);
                }
            }
        }
    }

    public double[] getBins() {
        return bins;
    }

    public void setPopulationFit(double populationFit) {
        this.populationFit = populationFit;
    }

    public double getPopulationFit() {
        return populationFit;
    }

    public String getPopulationName() {
        return populationName;
    }

    public static class ConfidenceInterval {

        private double start;
        private double stop;
        private double contamination;
        private double maxLogLikelihood;
        final double[] newBins;

        public ConfidenceInterval(double bins[], double intervalArea) {
            // make a copy of the bins in non-log space
            int maxIndex = MathUtils.maxElementIndex(bins);
            maxLogLikelihood = bins[maxIndex];

            newBins = MathUtils.applyToArray(bins, x -> Math.pow(10,x - bins[maxIndex]));
            final double total = MathUtils.sum(newBins)
            MathUtils.normalizeFromRealSpaceInPlace(newBins);

            double areaUnderCurve = 0;
            int leftIndex = maxIndex;
            int rightIndex = maxIndex;
            while (areaUnderCurve < 0.95) {

                // if the "left" bin is bigger, and can be moved, move it
                if (newBins[leftIndex] >= newBins[rightIndex] && leftIndex > 0) {
                    leftIndex--;
                } else {
                    // otherwise move the right bin if possible
                    if (rightIndex < bins.length - 1) {
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
                    areaUnderCurve += newBins[x];
            }
            start = (bins.length - rightIndex) * (100.0/bins.length);
            stop = (bins.length - leftIndex) * (100.0/bins.length);
            contamination = (bins.length - maxIndex) * (100.0/bins.length);
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

        public String toString() {
            return String.format("%f[%f-%f] log likelihood = %f", contamination, start, stop, maxLogLikelihood);
        }
    }
}
