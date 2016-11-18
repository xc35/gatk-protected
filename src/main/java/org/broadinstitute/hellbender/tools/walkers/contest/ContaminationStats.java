package org.broadinstitute.hellbender.tools.walkers.contest;

import org.broadinstitute.hellbender.utils.Nucleotide;

/**
 * a class that tracks our contamination stats; both the estimate of contamination, as well as the number of sites and other
 * run-specific data
 */
public final class ContaminationStats {
    private static final int ALLELE_COUNT = 4;
    private int numberOfSites = 0;
    private double sumOfAlleleFrequency = 0.0;
    private long basesFor = 0L;
    private long basesAgainst = 0L;
    private long basesOther = 0L;
    private ContaminationEstimate contaminationEstimate;
    private final Nucleotide.Counter alleleBreakdown;

    public ContaminationStats(int numberOfSites, double sumOfAlleleFrequency, long basesFor, long basesAgainst, long basesOther, Nucleotide.Counter alleleBreakdown, ContaminationEstimate estimate) {
        this.numberOfSites = numberOfSites;
        this.sumOfAlleleFrequency = sumOfAlleleFrequency;
        this.basesFor = basesFor;
        this.basesAgainst = basesAgainst;
        this.contaminationEstimate = estimate;
        this.alleleBreakdown = alleleBreakdown;
    }

    public long getBasesMatching() {
        return basesFor;
    }

    public long getBasesMismatching() {
        return basesAgainst;
    }

    public ContaminationEstimate getContamination() {
        return contaminationEstimate;
    }

    public void add(ContaminationStats other) {
        if (other == null) return;
        this.numberOfSites          += other.numberOfSites;
        this.sumOfAlleleFrequency   += other.sumOfAlleleFrequency;
        this.basesOther             += other.basesOther;
        this.basesFor               += other.basesFor;
        this.basesAgainst           += other.basesAgainst;
        this.alleleBreakdown.increment(other.alleleBreakdown);
        for (int i = 0; i < this.contaminationEstimate.getBins().length; i++) {
            this.contaminationEstimate.getBins()[i] += other.contaminationEstimate.getBins()[i];
        }
        this.contaminationEstimate.setPopulationFit(this.contaminationEstimate.getPopulationFit() +other.contaminationEstimate.getPopulationFit());
    }
}