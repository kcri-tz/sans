#include "graph.h"

/**
 * This is the size of the top list.
 */
uint64_t graph::t;

/**
 * This is a hash table mapping k-mers to colors [O(1)].
 */
unordered_map<kmer_t, color_t> graph::kmer_table;

/**
 * This is a hash table mapping colors to weights [O(1)].
 */
unordered_map<color_t, array<uint32_t,2>> graph::color_table;

/**
 * This is an ordered tree collecting the splits [O(log n)].
 */
multimap<double, color_t, greater<>> graph::split_list;

/**
 * This function initializes the top list size.
 *
 * @param t top list size
 */
void graph::init(uint64_t& top_size) {
    t = top_size;
}

/**
 * This function extracts k-mers from a sequence and adds them to the hash table.
 *
 * @param str dna sequence
 * @param color color flag
 */
void graph::add_kmers(string& str, uint64_t& color) {

    kmer_t kmer;    // create a new empty bit sequence for the k-mer
    kmer_t rcmer;    // create a bit sequence for the reverse complement

    next_kmer:
    if (str.length() < kmer::k) {
        return;    // not enough characters for a k-mer, abort
    }
    uint64_t pos = 0;    // current position in the string, from 0 to length

    for (; pos < kmer::k; ++pos) {    // collect the first k bases from the string
        if (str[pos] != 'A' && str[pos] != 'C' && str[pos] != 'G' && str[pos] != 'T') {
            str = str.substr(pos+1, string::npos);
            goto next_kmer;    // unknown base, start a new k-mer from the beginning
        }
        kmer::shift_right(kmer, str[pos]);    // shift each base into the bit sequence
    }
    rcmer = kmer;
    kmer::reverse_complement(rcmer, true);    // invert the k-mer, if necessary
    color::set(kmer_table[rcmer], color);    // update the k-mer with the current color

    for (; pos < str.length(); ++pos) {    // collect the remaining bases from the string
        if (str[pos] != 'A' && str[pos] != 'C' && str[pos] != 'G' && str[pos] != 'T') {
            str = str.substr(pos+1, string::npos);
            goto next_kmer;    // unknown base, start a new k-mer from the beginning
        }
        kmer::shift_right(kmer, str[pos]);    // shift each base into the bit sequence

        rcmer = kmer;
        kmer::reverse_complement(rcmer, true);    // invert the k-mer, if necessary
        color::set(kmer_table[rcmer], color);    // update the k-mer with the current color
    }
}

/**
 * This function iterates over the hash table and calculates the split weights.
 *
 * @param mean weight function
 */
void graph::add_weights(double mean(uint32_t&, uint32_t&)) {

    double min_value = numeric_limits<double>::min();    // current min. weight in the top list (>0)
    for (auto it = kmer_table.begin(); it != kmer_table.end(); ++it) {    // iterate over k-mer hash table
        color_t& color = it->second;    // get the color set for each k-mer
        bool pos = color::complement(color, true);    // invert the color set, if necessary
        array<uint32_t,2>& weight = color_table[color];    // get the weight and inverse weight for the color set

        double old_value = mean(weight[0], weight[1]);    // calculate the old mean value
        if (old_value >= min_value) {    // if it is greater than the min. value, find it in the top list
            auto range = split_list.equal_range(old_value);    // get all color sets with the given weight
            for (auto it = range.first; it != range.second; ++it) {
                if (it->second == color) {    // iterate over the color sets to find the correct one
                    split_list.erase(it);    // erase the entry with the old weight
                    break;
                }
            }
        }
        weight[pos]++;    // update the weight or the inverse weight of the current color set

        double new_value = mean(weight[0], weight[1]);    // calculate the new mean value
        if (new_value >= min_value) {    // if it is greater than the min. value, add it to the top list
            split_list.emplace(new_value, color);    // insert it at the correct position ordered by weight
            if (t != 0 && split_list.size() > t) {
                split_list.erase(--split_list.end());    // if the top list exceeds its limit, erase the last entry
                min_value = split_list.rbegin()->first;    // update the min. value for the next iteration
            }
        }
    }
}

/**
 * This function filters a maximum weight 1-tree compatible subset.
 */
void graph::filter_tree1() {
    auto tree1 = vector<color_t>();    // create a set for compatible splits
    for (auto it = split_list.begin(); it != split_list.end(); ) {
        if (test_set(it->second, tree1)) {
            tree1.emplace_back(it->second);
            ++it;    // if compatible, add the new split to the set
        } else {
            it = split_list.erase(it);    // otherwise, remove split
        }
    }
}

/**
 * This function filters a maximum weight 2-tree compatible subset.
 */
void graph::filter_tree2() {
    auto tree1 = vector<color_t>();    // create a set for compatible splits
    auto tree2 = vector<color_t>();    // create a set for incompatible splits
    for (auto it = split_list.begin(); it != split_list.end(); ) {
        if (test_set(it->second, tree1)) {
            tree1.emplace_back(it->second);
            ++it;    // if compatible, add the new split to the set
        } else if (test_set(it->second, tree2)) {
            tree2.emplace_back(it->second);
            ++it;    // if 2-tree compatible, add it to the backup set
        } else {
            it = split_list.erase(it);    // otherwise, remove split
        }
    }
}

/**
 * This function does not filter the splits in the output list.
 */
void graph::filter_none() {
    // split_list is not modified
}

/**
 * This function tests if a split is compatible with an existing set of splits.
 *
 * @param color new split
 * @param color_set set of splits
 * @return true, if compatible
 */
bool graph::test_set(color_t& color, vector<color_t>& color_set) {
    for (auto elem : color_set) {
        if (!color::is_compatible(color, elem)) {
            return false;    // compare to each split in the set
        }
    }
    return true;
}

/**
 * This function prints the splits ordered by weight to an output file stream.
 *
 * @param out output stream
 * @param names file names
 */
void graph::output_splits(ostream& out, vector<string>& names) {

    uint64_t pos = 0;
    for (auto& split : split_list) {
        out << split.first;    // weight of the split
        for (auto& name : names) {
            if (color::test(split.second, pos)) {
                out << "\t" << name;    // name of the file
            }
            split.second >>= 01u;
        }
        out << endl;
    }
}
