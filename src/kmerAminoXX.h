#include <iostream>
#include <bitset>
#include "util.h"

using namespace std;

/**
 * This is the max. k-mer length defined as a preprocessor directive.
 */
#ifndef maxK
#define maxK 32
#endif

/**
 * This class contains functions for working with k-mers containing amino acids.
 */
class kmerAminoXX {

private:

    /**
     * This is a bit-mask to erase all bits that exceed the k-mer length.
     */
    static bitset<5*maxK> mask;

public:

    /**
     * This is the length of a k-mer.
     */
    static uint64_t k;

    /**
    * The mod to use for binning
    */
    static uint64_t table_count;

    /**
    * The current binning carry
    */
    static uint64_t bin;

    /**
     * This function initializes the k-mer length and bit-mask.
     *
     * @param kmer_length k-mer length
     * @param bins the number of hash tables
     */
    static void init(uint64_t& kmer_length, uint64_t& bins);

    /**
     * This function shifts a k-mer adding a new character to the left.
     *
     * @param kmer bit sequence
     * @param c left character
     * @return right character
     */
    static char shift_left(bitset<5*maxK>& kmer, char& c);

    /**
     * This function shifts a k-mer adding a new character to the right.
     *
     * @param kmer bit sequence
     * @param c right character
     * @return left character
     */
    static char shift_right(bitset<5*maxK>& kmer, char& c);

    /**
     * This function constructs the reverse complement of a given k-mer.
     *
     * @param kmer bit sequence
     * @param minimize only invert, if smaller
     * @return 1 if inverted, 0 otherwise
     */
    static bool reverse_complement(bitset<5*maxK>& kmer, bool minimize);

    /**
     * This function computes the modulo on the Amino bitset
     *
     */
    static uint64_t bit_mod(const bitset<5*maxK>& kmer, uint64_t& module);


protected:

};
