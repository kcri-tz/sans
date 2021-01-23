#ifndef SANS_CLEANLINESS_H
#define SANS_CLEANLINESS_H

#include "graph.h"

/**
 * This struct saves unfiltered and filtered split values to
 * compare them after filtering so we don't need to save the entire unfiltered list
 */
struct filter_weight_stats {
    double weight;
    color_t color;

    /**
     * checks, if the given other struct is the same as this one
     * @param other
     * @return true, if weight and color are the same
     */
    bool is(filter_weight_stats other) {
        return other.weight == weight && other.color == color;
    }
};

/**
 *  Cleanliness provides functions for saving, calculating and reporting the cleanliness ratio
 *  after the filter has been applied.
 */
class cleanliness {

private:
    /**
     * stores the index of the last split after the filtering process
     */
    uint64_t splitCountAfter = 0;
    /**
    * stores the weight sum after the filtering process
    */
    double splitWeightCountAfter = 0;

    /**
     * stores the index of the last split before the filtering process
     */
    uint64_t splitCountBefore = 0;
    /**
    * stores the weight sum before the filtering process
    */
    double splitWeightCountBefore = 0;

    /**
    * stores the smallest weight stat for determine the last counted split
    */
    struct filter_weight_stats smallestWeight = {};

    /**
    * stores the skipped weight stats before the filtering
    */
    vector<filter_weight_stats> weightsBefore;

public:
    /**
     * Currently this method does nothing and is only for future expansion purpose
     */
    void init();

    /**
     * adds a Weight State before the filtering
     */
    void addWeightStateBefore(double weight, color_t color);

    /**
     * sets a new smallest weight after processing the filtered split
     */
    void setSmallestWeight(double weight, color_t color);

    /**
     * increments the after counters with the given weight
     */
    void incrementWeightAfterCounter(double weight);

    /**
     * calculates the Before counters in dependence of the given after counters
     */
    void calculateWeightBeforeCounter();

    /**
     * sets the size of the splits after filtering
     */
    void setFilteredCount(uint64_t splitCountAfter);

    /**
     * prints the state of the cleanliness object
     */
    void reportCleanliness();
};


#endif //SANS_CLEANLINESS_H
