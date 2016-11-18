package org.broadinstitute.hellbender.tools.walkers.contest;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.*;

/**
 * our contamination results object; this object aggregates the results of the contamination run over lanes, samples,
 * or whatever other divisor we've used on the read data
 */
public final class ContaminationResults {
    protected static final Logger logger = LogManager.getLogger(ContaminationResults.class);


    // a map of our contamination targets and their stats
    // key: aggregation entity ("BAM", sample name, or lane name)
    // value: ContaminationStats (whcih
    private Map<String, ContaminationEstimate> stats = new HashMap<>();

    final Map<String, List<ContaminationEstimate>> storedData = new HashMap<>();

    /**
     * add to the stats
     *
     * @param newAggregationStats a mapping of the stat name to their statistics collected
     */
    public void add(final Map<String, ContaminationEstimate> newAggregationStats) {

        // for each aggregation level
        for (final String aggregationKey : newAggregationStats.keySet()) {
            final ContaminationEstimate newStats = newAggregationStats.get(aggregationKey);

            // a new way of doing this... store all the data until the end...
            if (!storedData.containsKey(aggregationKey)) {
                storedData.put(aggregationKey, new ArrayList<>());
            }

            storedData.get(aggregationKey).add(newStats);

            // merge the sets
            if (stats.containsKey(aggregationKey)) {
                stats.get(aggregationKey).add(newStats);
            } else {
                stats.put(aggregationKey, newStats);
            }
        }
    }


    public void outputReport(final double precision) {

        //TODO: logger.info isn't correct -- use a TableUtils method
        logger.info("name\tpopulation\tpopulation_fit\tcontamination\tconfidence_interval_95_width\tconfidence_interval_95_low\tconfidence_interval_95_high\tsites");

        //TODO: very confused why stats = entry.getValue isn't used
        // TODO: and why have both storedData and stats???
        for (final Map.Entry<String, ContaminationEstimate> entry : stats.entrySet()) {
            final ContaminationEstimate stats = entry.getValue();
            final String aggregationLevel = entry.getKey();

            final List<ContaminationEstimate> data = storedData.get(aggregationLevel);
            final String pm = "%3." + Math.round(Math.log10(1/precision)) +"f";

            final int bins = data.iterator().next().getBins().length;

            //TODO: sort before output?

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
            logger.info(String.format("%s\t%s\t"+pm+"\t"+pm+"\t"+pm+"\t"+pm+"\t"+"%d",
                    aggregationLevel,
                    "n/a",
                    newInterval.getContamination(),
                    (newInterval.getStop() - newInterval.getStart()),
                    newInterval.getStart(),
                    newInterval.getStop(),
                    data.size()));
        }
    }
}