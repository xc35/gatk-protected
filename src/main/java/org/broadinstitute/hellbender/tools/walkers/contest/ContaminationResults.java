package org.broadinstitute.hellbender.tools.walkers.contest;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

/**
 * our contamination results object; this object aggregates the results of the contamination run over lanes, samples,
 * or whatever other divisor we've used on the read data
 */
public final class ContaminationResults {
    protected static final Logger logger = LogManager.getLogger(ContaminationResults.class);

    public static class ContaminationData {
        private final double[] bins;
        private final double p;

        public double[] getBins() {
            return bins;
        }

        public double getP() {
            return p;
        }

        public ContaminationData(final long basesMatching, final long basesMismatching, final double[] bins) {
            this.bins = bins;

            final BetaDistribution dist = new BetaDistribution(basesMismatching + 1, basesMatching + 1);
            this.p = 1.0d - dist.cumulativeProbability(0.5d);
        }
    }

    // a map of our contamination targets and their stats
    // key: aggregation entity ("BAM", sample name, or lane name)
    // value: ContaminationStats (whcih
    private Map<String,Map<String, ContaminationStats>> stats = new HashMap<String,Map<String, ContaminationStats>>();

    final Map<String, Map<String, List<ContaminationData>>> storedData = new HashMap<String, Map<String, List<ContaminationData>>>();

    /**
     * add to the stats
     *
     * @param newAggregationStats a mapping of the stat name to their statistics collected
     */
    public void add(final Map<String, Map<String, ContaminationStats>> newAggregationStats) {

        // for each aggregation level
        for (final String aggregationKey : newAggregationStats.keySet()) {
            final Map<String, ContaminationStats> contaminationStatsByPopulation = newAggregationStats.get(aggregationKey);

            // a new way of doing this... store all the data until the end...
            if (!storedData.containsKey(aggregationKey)) {
                storedData.put(aggregationKey, new HashMap<>());
            }
            for (final String pop : contaminationStatsByPopulation.keySet()) {
                final ContaminationStats newStats = contaminationStatsByPopulation.get(pop);

                // if it exists... just merge it
                if (!storedData.get(aggregationKey).containsKey(pop)) {
                    storedData.get(aggregationKey).put(pop, new ArrayList<>());
                }

                final double[] newData = new double[newStats.getContamination().getBins().length];
                System.arraycopy(newStats.getContamination().getBins(),0,newData,0,newStats.getContamination().getBins().length);
                storedData.get(aggregationKey).get(pop).add(new ContaminationData(newStats.getBasesMatching(), newStats.getBasesMismatching(), newData));
            }

            // merge the sets
            if (stats.containsKey(aggregationKey)) {

                // and for each population
                for (final String pop : contaminationStatsByPopulation.keySet()) {
                    final ContaminationStats newStats = contaminationStatsByPopulation.get(pop);

                    // if it exists... just merge it
                    if (stats.get(aggregationKey).containsKey(pop)) {
                        stats.get(aggregationKey).get(pop).add(newStats);
                    } else {
                        stats.get(aggregationKey).put(pop, newStats);
                    }
                }
            } else {
                stats.put(aggregationKey, contaminationStatsByPopulation);
            }
        }
    }

    /**
     * output the contamination data, and return the contamination data
     */
    public void outputReport(final double precision, final double betaThreshold) {

        //TODO: logger.info isn't correct -- use a TableUtils method
        logger.info("name\tpopulation\tpopulation_fit\tcontamination\tconfidence_interval_95_width\tconfidence_interval_95_low\tconfidence_interval_95_high\tsites");

        for (final Map.Entry<String,Map<String, ContaminationStats>> entry : stats.entrySet()) {
            for (final ContaminationStats stats : entry.getValue().values()) {
                final String aggregationLevel = entry.getKey();
                final String population = stats.getContamination().getPopulationName();

                final List<ContaminationData> newStats = storedData.get(aggregationLevel).get(population);
                final String pm = "%3." + Math.round(Math.log10(1/precision)) +"f";

                final int bins = newStats.iterator().next().getBins().length;

                // sort the collection
                Collections.sort(newStats);

                final List<ContaminationData> data = newStats.stream().filter(d -> d.getP() < betaThreshold).collect(Collectors.toList());


                final double[][] matrix = new double[bins][data.size()];

                for (int i = 0; i<bins; i++) {
                    for (int j=0; j<data.size(); j++) {
                        matrix[i][j] = data.get(j).getBins()[i];
                    }
                }

                // now perform the sum
                final double[] output = new IndexRange(0, bins).mapToDouble(n -> MathUtils.sum(matrix[n]));

                // get the confidence interval, at the set width
                final ContaminationEstimate.ConfidenceInterval newInterval = new ContaminationEstimate.ConfidenceInterval(output);

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

    public void writeCurves(final PrintStream out) {
        boolean outputBins = false;
        for (final Map.Entry<String, Map<String, ContaminationStats>> entry : stats.entrySet()) {
            for (final ContaminationStats stats : entry.getValue().values()) {
                if (!outputBins) {
                    final String[] bins = new String[stats.getContamination().getBins().length];
                    for (int index = 0; index < stats.getContamination().getBins().length; index++) {
                        bins[index] = String.valueOf(100.0 * (1 - (double) index / stats.getContamination().getBins().length));
                    }
                    outputBins = true;
                    out.print("name,pop,");
                    out.println(Utils.join(",",(Object)bins));
                }
                final String[] bins = new String[stats.getContamination().getBins().length];
                int index = 0;
                for (final double value : stats.getContamination().getBins()) {
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

    public void setStats(final Map<String, Map<String,ContaminationStats>> stats) {
        this.stats = stats;
    }
}