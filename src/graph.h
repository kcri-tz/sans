#include <iostream>

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <limits>

using namespace std;

#include "kmer32.h"
#include "kmerXX.h"

#if maxK > 32 // store k-mers in a bitset, allows larger k-mers
    typedef kmerXX kmer;
    typedef bitset<2*maxK> kmer_t;
#else // store k-mer bits in an integer, optimizes performance
    typedef kmer32 kmer;
    typedef uint64_t kmer_t;
#endif

#include "color64.h"
#include "colorXX.h"

#if maxN > 64 // store colors in a bitset, allows more input files
    typedef colorXX color;
    typedef bitset<maxN> color_t;
#else // store color bits in an integer, optimizes performance
    typedef color64 color;
    typedef uint64_t color_t;
#endif

/**
 * This class manages the k-mer/color hash tables and split list.
 */
class graph {

private:

    /**
     * This is a hash table mapping k-mers to colors [O(1)].
     */
    static unordered_map<kmer_t, color_t> kmer_table;

    /**
     * This is a hash table mapping colors to weights [O(1)].
     */
    static unordered_map<color_t, array<uint32_t,2>> color_table;

public:

    /**
     * This is an ordered tree collecting the splits [O(log n)].
     */
    static multimap<double, color_t, greater<>> split_list;

    /**
     * This is the size of the top list.
     */
    static uint64_t t;

    /**
     * This function initializes the top list size.
     *
     * @param t top list size
     */
    static void init(uint64_t& top_size);

    /**
     * This function extracts k-mers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     */
    static void add_kmers(string& str, uint64_t& color);

    /**
     * This function extracts k-mers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param max_iupac allowed number of ambiguous k-mers per position
     */
    static void add_kmers(string& str, uint64_t& color, uint64_t& max_iupac);

    /**
     * This function iterates over the hash table and calculates the split weights.
     *
     * @param mean weight function
     */
    static void add_weights(double mean(uint32_t&, uint32_t&));

    /**
     * This function adds a single split (weight and colors) to the output list.
     *
     * @param weight split weight
     * @param color split colors
     */
    static void add_split(double& weight, color_t& color);

    /**
     * This function filters a greedy maximum weight tree compatible subset.
     *
     * @param verbose print progress
     */
    static void filter_strict(bool& verbose);

    /**
     * This function filters a greedy maximum weight weakly compatible subset.
     *
     * @param verbose print progress
     */
    static void filter_weakly(bool& verbose);

    /**
     * This function filters a greedy maximum weight n-tree compatible subset.
     *
     * @param n number of trees
     * @param verbose print progress
     */
    static void filter_n_tree(uint64_t n, bool& verbose);

protected:

    /**
     * This function tests if a split is compatible with an existing set of splits.
     *
     * @param color new split
     * @param color_set set of splits
     * @return true, if compatible
     */
    static bool test_strict(color_t& color, vector<color_t>& color_set);

    /**
     * This function tests if a split is weakly compatible with an existing set of splits.
     *
     * @param color new split
     * @param color_set set of splits
     * @return true, if weakly compatible
     */
    static bool test_weakly(color_t& color, vector<color_t>& color_set);

    /**
     * This function calculates the multiplicity of iupac k-mers.
     *
     * @param product overall multiplicity
     * @param factors per base multiplicity
     * @param input iupac character
     */
    static void iupac_calc(long double& product, vector<uint8_t>& factors, char& input);

    /**
     * This function shifts a base into a set of ambiguous iupac k-mers.
     *
     * @param prev set of k-mers
     * @param next set of k-mers
     * @param input iupac character
     */
    static void iupac_shift(unordered_set<kmer_t>& prev, unordered_set<kmer_t>& next, char& input);

};
