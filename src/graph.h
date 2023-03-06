#include <iostream>
#include <limits>
#include <algorithm>
#include <functional>
#include <utility>
#include <vector>
#include <thread>

#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "tsl/sparse_map.h"
#include "tsl/sparse_set.h"

using namespace std;

template <typename K, typename V>
    // using hash_map = unordered_map<K,V>;
    using hash_map = tsl::sparse_pg_map<K,V>;
template <typename T>
    // using hash_set = unordered_set<T>;
    using hash_set = tsl::sparse_pg_set<T>;

#include "kmer32.h"
#include "kmerXX.h"
#include "kmerAmino12.h"
#include "kmerAminoXX.h"

#if maxK > 32 // store k-mers in a bitset, allows larger k-mers
    typedef kmerXX kmer;
    typedef bitset<2*maxK> kmer_t;
#else // store k-mer bits in an integer, optimizes performance
    typedef kmer32 kmer;
    typedef uint64_t kmer_t;
#endif

#if maxK > 12 // store k-mers in a bitset, allows larger k-mers
    typedef kmerAminoXX kmerAmino;
    typedef bitset<5*maxK> kmerAmino_t;
#else // store k-mer bits in an integer, optimizes performance
    typedef kmerAmino12 kmerAmino;
    typedef uint64_t kmerAmino_t;
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

#include "spinlockMutex.h"


/**
 * A tree structure that is needed for generating a NEWICK string.
 */
struct node {
    color_t taxa;
    double weight;
    vector<node*> subsets;
};


/**
 * This class manages the k-mer/color hash tables and split list.
 */
class graph {

private:

    /**
     * This is the size of the top list.
     */
    static uint64_t t;

    /**
     * This is the min. coverage threshold for k-mers.
     */
    static uint64_t quality;

    static bool isAmino;
    
    /**
     * This int indicates the number of tables to use for hashing
     */
    static uint64_t table_count;
    
    /**
     * This vector holds the carries of 2**i % table_count for fast distribution of bitset represented kmers
     */
    static vector<uint64_t> period;
    static uint64_t first_mod_correction;
    static uint64_t second_mod_correction;
    /**
     * This is a vector of hash tables mapping k-mers to colors [O(1)].
     */
    static vector<hash_map<kmer_t, color_t>> kmer_table;

    /**
     * This is a vector of spinlocks protecting the hash tables.
     */
    static vector<std::mutex> lock;

    /**
     * This is a hash table mapping k-mers to colors [O(1)].
     */
    static vector<hash_map<kmerAmino_t, color_t>> kmer_tableAmino;

    /**
     * This is a hash table mapping colors to weights [O(1)].
     */
    static hash_map<color_t, array<uint32_t,2>> color_table;

    /**
     * This is a hash set used to filter k-mers for coverage (q > 1).
     */
    static vector<hash_set<kmer_t>> quality_set;
    static vector<hash_set<kmerAmino_t>> quality_setAmino;

    /**
     * This is a hash map used to filter k-mers for coverage (q > 2).
     */
    static vector<hash_map<kmer_t, uint64_t>> quality_map;
    static vector<hash_map<kmerAmino_t, uint64_t>> quality_mapAmino;

public:
     
    // [Temporary: Test]
    static void showTableSizes();

    /**
     * This is an ordered tree collecting the splits [O(log n)].
     */
    static multimap<double, color_t, greater<>> split_list;

    /**
    * These are the allowed chars.
    */
    static vector<char> allowedChars;

    /**
     * This function initializes the top list size, coverage threshold, and allowed characters.
     *
     * @param top list size
     * @param isAmino use amino processing
     * @param quality coverage threshold
     * @param bins hash_tables to use for parallel processing
     * @param thread_count the number of threads used for processing
     */
    static void init(uint64_t& top_size, bool isAmino, uint64_t& quality, uint64_t& bins, uint64_t& thread_count);

    /**
     * This function shift updates the bin of a kmer
    */
    static uint64_t shift_update_bin(uint64_t& bin, char& left, char& right);
    
    /**
     * This method shift updates the reverse complement bin for a kmer
    */
   static uint64_t shift_update_rc_bin(uint64_t& rc_bin, char& c_left, char& c_right);

    /**
    * This function shift updates a bin for an amino kmer
    */
    static uint64_t shift_update_amino_bin(uint64_t bin, kmerAmino_t& kmer, char& c_left, char& c_right);

    /**
     *  This method computes the bin of a given kmer(slower than shift update)
     * @param kmer The target kmer
     * @return uint64_t The bin
     */
    #if (maxK <= 32)
    static uint64_t compute_bin(const kmer_t& kmer);
    #else
    static uint64_t compute_bin(const bitset<2*maxK>& kmer);
    #endif

    /**
     *  This function computes the bin of a given amino kmer(slower than shift update)
     * @param kmer The target kmer
     * @return uint64_t The bin
     */
    #if (maxK <= 12)
        static uint64_t compute_amino_bin(const kmerAmino_t kmer);
    #else
        static uint64_t compute_amino_bin(const bitset<5*maxK>& kmer);
    #endif


    /**
     * This function computes the target hash table index for a given k-mer
     * @param k-mer The k-mer to compute the index for
     * @param reversed The bool implying if the k-mer was reversed of not 
     */
    static uint64_t get_table_index(const kmer_t& kmer, bool reversed);

    /**
     * This function computes the target hash table index for a given amino k-mer
     * @param kme The k-mer to compute the index for
     *
     */
    static uint64_t get_amino_table_index(const kmerAmino_t& kmer);

    /**
     * This function hashes a base k-mer and stores it in the corresponding hash table
     *  @param t_id The thread_id of the current thread
     *  @param kmer The k-mer to store
     *  @param color The color to store 
     *  @param reversed The bool implying if the k-mer was reversed or not
     */
    static void hash_kmer(uint64_t bin, const kmer_t& kmer, const uint64_t& color);


    /**
     * This function hashes a base k-mer and stores it in the corresponding hash table (sequential version)
     *  @param kmer The k-mer to store
     *  @param color The color to store 
     *  @param reversed The bool implying if the k-mer was reversed or not
     */
    static void hash_kmer(const kmer_t& kmer, const uint64_t& color);



    /**
     * This function hashes an amino k-mer and stores it in the correstponding hash table
     *  @param t_id The id of the current thread
     *  @param kmer The kmer to store
     *  @param color The color to store 
     */
    static void hash_kmer_amino(uint64_t bin, const kmerAmino_t& kmer, const uint64_t& color);

    /**
     * This function hashes an amino k-mer and stores it in the correstponding hash table (sequential version)
     *  @param kmer The kmer to store
     *  @param color The color to store 
     */
    static void hash_kmer_amino(const kmerAmino_t& kmer, const uint64_t& color);


    /**
     * This function searches the bit-wise corresponding hash table for the given kmer
     * @param kmer The kmer to search
     * @return True if the kmer is stored.
     */
    static bool search_kmer(const kmer_t& kmer);

    /**
     * This function searches the bit-wise corresponding hash table for the given amnio kmer
     * @param kmer The kmer to search
     * @return Ture if the kmer is stored
     */
    static bool search_kmer_amino(const kmerAmino_t& kmer);


    /**
     * This function returns the stored colores of the given kmer
     * @param kmer The target kmer
     * @return color_t The stored colores
     */
    static color_t get_color(const kmer_t& kmer, bool reversed);

    /**
     * This function returns the stored color of the given kmer
     * @param kmer The target kmer
     * @return color_t The stored color
     */
    static color_t get_color_amino(const kmerAmino_t& kmer);


    /**
     * This function removes the kmer entry from the hash map
     */
    static void remove_kmer(const kmer_t& kmer, bool reversed);

    /**
     * This function removes the amino kmer entry from the corresponding hash table
     */
    static void remove_kmer_amino(const kmerAmino_t& kmer);

    /**
     * This function extracts k-mers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param reverse merge complements
     */
    static void add_kmers(uint64_t& T, string& str, uint64_t& color, bool& reverse);

    /**
     * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param reverse merge complements
     * @param m number of k-mers to minimize
     */
    static void add_minimizers(uint64_t& T, string& str, uint64_t& color, bool& reverse, uint64_t& m);

    /**
     * This function extracts k-mers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param reverse merge complements
     * @param max_iupac allowed number of ambiguous k-mers per position
     */
    static void add_kmers(uint64_t& T, string& str, uint64_t& color, bool& reverse, uint64_t& max_iupac);

    /**
     * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
     *
     * @param str dna sequence
     * @param color color flag
     * @param reverse merge complements
     * @param m number of k-mers to minimize
     * @param max_iupac allowed number of ambiguous k-mers per position
     */
    static void add_minimizers(uint64_t& T, string& str, uint64_t& color, bool& reverse, uint64_t& m, uint64_t& max_iupac);

    /**
     * This function calculates the split weight for a single entry of the hash table
     * 
     * @param mean weight function
     * @param verbose print progress
     * @param min_value the minimal weight represented in the top list
     * @return the new minimal weight represented in the top list
     */
    static double add_weight(color_t& color, double mean(uint32_t&, uint32_t&), double min_value, bool pos);

    /**
     * This function iterates over the hash table and calculates the split weights.
     *
     * @param mean weight function
     * @param verbose print progress
     * @param min_value the minimal weight currently represented in the top list
     */
    static void add_weights(double mean(uint32_t&, uint32_t&), double min_value, bool& verbose);

    /**
     * This function adds a single split (weight and colors) to the output list.
     *
     * @param weight split weight
     * @param color split colors
     */
    static void add_split(double& weight, color_t& color);

    /**
     * This funtion adds a sigle split from a cdbg to the output list.
     * 
     * @param mean mean function
     * @param kmer_seq kmer
     * @param kmer_color split colors
     * @param min_value the minimal weight currently represented in the top list
     * @return The minimal weight currently represented in the top list
     */
     static double add_cdbg_colored_kmer(double mean(uint32_t&, uint32_t&), string kmer_seq, color_t& kmer_color, double min_value);       

    /**
     * This function clears color-related temporary files.
     */
    static void clear_thread(uint64_t& T);

    /**
     * This function filters a greedy maximum weight tree compatible subset.
     *
     * @param verbose print progress
     * @return the new minimal weight represented in the top list
     */
    static void filter_strict(bool& verbose);

    /**
     * This function filters a greedy maximum weight tree compatible subset and returns a newick string.
     *
     * @param map function that maps an integer to the original id, or null
     * @param verbose print progress
     */
    static string filter_strict(std::function<string(const uint64_t&)> map, bool& verbose);

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

    /**
     * This function filters a greedy maximum weight n-tree compatible subset and returns a string with all trees in newick format.
     *
     * @param n number of trees
     * @param map function that maps an integer to the original id, or null
     * @param verbose print progress
     */
    static string filter_n_tree(uint64_t n, std::function<string(const uint64_t&)> map, bool& verbose);

protected:

    /**
     * This function qualifies a k-mer and places it into the hash table.
     *
     * @param kmer bit sequence
     * @param color color flag
     */
    static function<void(uint64_t& T, uint64_t& bin, const kmer_t&, uint64_t&)> emplace_kmer;

    /**
     * This function qualifies a k-mer and places it into the hash table.
     *
     * @param kmer bit sequence
     * @param color color flag
     */
    static function<void(uint64_t& T, uint64_t& bin, const kmerAmino_t&, uint64_t&)> emplace_kmer_amino;

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
    static void iupac_shift(hash_set<kmer_t>& prev, hash_set<kmer_t>& next, char& input);

    /**
   * This function shifts a base into a set of ambiguous iupac k-mers.
   *
   * @param prev set of k-mers
   * @param next set of k-mers
   * @param input iupac character
   */
    static void iupac_shift_amino(hash_set<kmerAmino_t>& prev, hash_set<kmerAmino_t>& next, char& input);

    /**
     * This function returns a tree structure (struct node) generated from the given list of color sets.
     *
     * @param color_set list of color sets
     * @return tree structure (struct node)
     */
    static node* build_tree(vector<color_t>& color_set);

    /**
     * This function recursively refines a given set/tree structure by a given split.
     *
     * @param current_set node of currently considered (sub-)set/tree structure
     * @param split color set to refine by
     * @return whether or not the given split is compatible with the set/tree structure
     */
    static bool refine_tree(node* current_set, color_t& split, color_t& allTaxa);

    /**
     * This function returns a newick string generated from the given tree structure (set).
     *
     * @param root root of the tree/set structure
     * @return newick string
     */
    static string print_tree(node* root, std::function<string(const uint64_t&)> map);

    /**
     * This function checks if the character at the given position is allowed.
     * @param pos position in str
     * @param str the current part of the sequence
     * @return true if allowed, false otherwise
     */
    static bool isAllowedChar(uint64_t pos, string &str);
};
