package org.broadinstitute.hellbender.tools.walkers.contest;

import htsjdk.samtools.util.Locatable;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.PrintStream;
import java.util.*;

/**
 * our contamination results object; this object aggregates the results of the contamination run over lanes, samples,
 * or whatever other divisor we've used on the read data
 */
public final class ContaminationResults {
    protected static final Logger logger = LogManager.getLogger(ContaminationResults.class);

    public static class ContaminationData implements Comparable<ContaminationData> {
        private Locatable site;
        private long basesMatching = 0l;
        private long basesMismatching = 0l;
        private double mismatchFraction = -1d;
        private double[] bins;
        private double p;

        public long getBasesMatching() {
            return basesMatching;
        }

        public long getBasesMismatching() {
            return basesMismatching;
        }

        public double getMismatchFraction() {
            return mismatchFraction;
        }

        public double[] getBins() {
            return bins;
        }

        public double getP() {
            return p;
        }

        public ContaminationData(Locatable site, long basesMatching, long basesMismatching, double[] bins) {
            this.site = site;
            this.basesMatching = basesMatching;
            this.basesMismatching = basesMismatching;
            this.bins = bins;
            long totalBases = this.basesMatching + this.basesMismatching;
            mismatchFraction = totalBases > 0 ? (double)this.basesMismatching / (double) totalBases : -1d;


            int a = (int) this.getBasesMismatching() + 1;
            int b = (int) this.getBasesMatching() + 1;
            BetaDistribution dist = new BetaDistribution(a,b);
            this.p = 1.0d - dist.cumulativeProbability(0.5d);
        }

        public int compareTo(ContaminationData other) {
            return -Double.compare(this.getP(), other.getP());
        }

        @Override
        public String toString() {
            return String.format("ContaminationData{site=%s, basesMatching=%d, basesMismatching=%d, mismatchFraction=%f}",site, basesMatching, basesMismatching, mismatchFraction);
        }
    }


    // what precision are we using in our calculations
    private final double precision;

    // a map of our contamination targets and their stats
    // key: aggregation entity ("META", sample name, or lane name)
    // value: ContaminationStats (whcih
    private Map<String,Map<String, ContaminationStats>> stats = new HashMap<String,Map<String, ContaminationStats>>();

    public ContaminationResults(double precision) {
        this.precision = precision;
    }


    Map<String, Map<String, List<ContaminationData>>> storedData = new HashMap<String, Map<String, List<ContaminationData>>>();

    /**
     * add to the stats
     *
     * @param newAggregationStats a mapping of the stat name to their statistics collected
     */
    public void add(final Map<String, Map<String, ContaminationStats>> newAggregationStats) {

        // for each aggregation level
        for (String aggregationKey : newAggregationStats.keySet()) {
            Map<String, ContaminationStats> populationContaminationStats = newAggregationStats.get(aggregationKey);

            // a new way of doing this... store all the data until the end...
            if (!storedData.containsKey(aggregationKey)) { storedData.put(aggregationKey, new HashMap<>()); }
            for (String pop : populationContaminationStats.keySet()) {
                ContaminationStats newStats = populationContaminationStats.get(pop);

                // if it exists... just merge it
                if (!storedData.get(aggregationKey).containsKey(pop)) {
                    storedData.get(aggregationKey).put(pop, new ArrayList<>());
                }

                double[] newData = new double[newStats.getContamination().getBins().length];
                System.arraycopy(newStats.getContamination().getBins(),0,newData,0,newStats.getContamination().getBins().length);
                storedData.get(aggregationKey).get(pop).add(new ContaminationData(newStats.getSite(), newStats.getBasesMatching(), newStats.getBasesMismatching(), newData));
            }

            // merge the sets
            if (stats.containsKey(aggregationKey)) {

                // and for each population
                for (String pop : populationContaminationStats.keySet()) {
                    ContaminationStats newStats = populationContaminationStats.get(pop);

                    // if it exists... just merge it
                    if (stats.get(aggregationKey).containsKey(pop)) {
                        stats.get(aggregationKey).get(pop).add(newStats);
                    } else {
                        stats.get(aggregationKey).put(pop, newStats);
                    }
                }
            } else {
                stats.put(aggregationKey, populationContaminationStats);
            }
        }
    }

    /**
     * output the contamination data, and return the contamination data
     * @return the contamination value
     */
    public void outputReport(double precision, double fractionToTrim, double trimInterval, double betaThreshold) {

        //TODO: logger.info isn't correct -- use a TableUtils method
        logger.info("name\tpopulation\tpopulation_fit\tcontamination\tconfidence_interval_95_width\tconfidence_interval_95_low\tconfidence_interval_95_high\tsites");

        for (Map.Entry<String,Map<String, ContaminationStats>> entry : stats.entrySet()) {
            for (ContaminationStats stats : entry.getValue().values()) {
                String aggregationLevel = entry.getKey();
                String population = stats.getContamination().getPopulationName();

                List<ContaminationData> newStats = storedData.get(aggregationLevel).get(population);
                String pm = "%3." + Math.round(Math.log10(1/precision)) +"f";

                int bins = newStats.iterator().next().getBins().length;
                int maxTrim = (int) Math.floor((double)(newStats.size()) * fractionToTrim);

                // sort the collection
                Collections.sort(newStats);

                List<ContaminationData> data = new ArrayList<ContaminationData>(newStats);

                // trim sites with > 95% p of being > 0.5 f (based on beta distribution)
                int trimmed = 0;
                for(Iterator<ContaminationData> i = data.iterator(); trimmed < maxTrim && i.hasNext();) {
                    ContaminationData x = i.next();
                    if (x.getP() >= betaThreshold) {
                        logger.info("Trimming " + x.toString() + " with p(f>=0.5) >= " + betaThreshold + " with a value of  " + x.getP());
                        i.remove();
                        trimmed++;
                    }
                }

                double[][] matrix = new double[bins][data.size()];

                for (int i = 0; i<bins; i++) {
                    for (int j=0; j<data.size(); j++) {
                        matrix[i][j] = data.get(j).getBins()[i];
                    }
                }

                // now perform the sum
                final double[] output = new IndexRange(0, bins).mapToDouble(n -> MathUtils.sum(matrix[n]));

                // get the confidence interval, at the set width
                ContaminationEstimate.ConfidenceInterval newInterval = new ContaminationEstimate.ConfidenceInterval(output, 0.95);

                //TODO: ditto about a Table method
                logger.info(String.format("%s\t%s\t%s\t"+pm+"\t"+pm+"\t"+pm+"\t"+pm+"\t"+"%d",
                        aggregationLevel,
                        population,
                        "n/a",
                        newInterval.getContamination(),
                        (newInterval.getStop() - newInterval.getStart()),
                        newInterval.getStart(),
                        newInterval.getStop(),
                        data.size()));
            }
        }
    }

    public void writeCurves(PrintStream out) {
        boolean outputBins = false;
        for (Map.Entry<String, Map<String, ContaminationStats>> entry : stats.entrySet()) {
            for (ContaminationStats stats : entry.getValue().values()) {
                if (!outputBins) {
                    String[] bins = new String[stats.getContamination().getBins().length];
                    for (int index = 0; index < stats.getContamination().getBins().length; index++) {
                        bins[index] = String.valueOf(100.0 * (1 - (double) index / stats.getContamination().getBins().length));
                    }
                    outputBins = true;
                    out.print("name,pop,");
                    out.println(Utils.join(",",(Object)bins));
                }
                String[] bins = new String[stats.getContamination().getBins().length];
                int index = 0;
                for (double value : stats.getContamination().getBins()) {
                    bins[index++] = String.valueOf(value);
                }
                out.print(entry.getKey()+",\""+stats.getContamination().getPopulationName()+"\",");
                out.println(Utils.join(",", (Object)bins));
            }
        }
    }

    public Map<String, Map<String, ContaminationStats>> getStats() {
        return Collections.unmodifiableMap(stats);
    }

    public void setStats(Map<String, Map<String,ContaminationStats>> stats) {
        this.stats = stats;
    }
}