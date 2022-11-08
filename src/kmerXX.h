#include <iostream>
#include <bitset>
#include <vector>

using namespace std;

/**
 * This is the max. k-mer length defined as a preprocessor directive.
 */
#ifndef maxK
#define maxK 32
#endif

/**
 * This class contains functions for working with k-mers (maxK > 32).
 */
class kmerXX {

private:

    /**
     * This is a bit-mask to erase all bits that exceed the k-mer length.
     */
    static bitset<2*maxK> mask;

public:

    /**
     * This is the length of a k-mer.
     */
    static uint64_t k;

    /**
     * The mod to use for binning
     */
    static uint64_t mod;

    /**
     * The hash table index of the kmer
     */
    static uint64_t bin;
    static uint64_t rbin;
    /**
     * This function initializes the k-mer length and bit-mask.
     *
     * @param kmer_length k-mer length
     * @param bins the number of hash tables to use for parallel hashing 
     */
    static void init(uint64_t& kmer_length, uint64_t& bins);


    /**
     * This function shifts a k-mer adding a new character to the left.
     *
     * @param kmer bit sequence
     * @param c left character
     * @return right character
     */
    static char shift_left(bitset<2*maxK>& kmer, char& c);

    /**
     * This function shifts a k-mer adding a new character to the right.
     *
     * @param kmer bit sequence
     * @param c right character
     * @return left character
     */
    static char shift_right(bitset<2*maxK>& kmer, char& c);

    /**
     * This function constructs the reverse complement of a given k-mer.
     *
     * @param kmer bit sequence
     * @param minimize only invert, if smaller
     * @return 1 if inverted, 0 otherwise
     */
    static bool reverse_complement(bitset<2*maxK>& kmer, bool minimize);

    /**
     * This function computes the modulo function on the bitset
     */
    static uint64_t bit_mod(const bitset<2*maxK>& kmer, uint64_t& module);

protected:

    /**
     * This function encodes a single character to two bits.
     *
     * @param c character
     * @return bit sequence
     */
    static uint64_t char_to_bits(char& c);

    /**
     * This function decodes two bits to a single character.
     *
     * @param b bit sequence
     * @return character
     */
    static char bits_to_char(uint64_t& b);

};
