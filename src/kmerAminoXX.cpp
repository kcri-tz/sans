#include "kmerAminoXX.h"

/**
 * This is the length of a k-mer.
 */
uint64_t kmerAminoXX::k;

/**
* The mod to use for binning
*/
uint64_t kmerAminoXX::table_count;

/**
 * The target hash table of the kmer
 */
uint64_t kmerAminoXX::bin;

/**
 * This vector stores the precomputed carries of 2^i % bin
 */
vector<uint64_t> kmerAminoXX::period;

/**
 * This is a bit-mask to erase all bits that exceed the k-mer length.
 */
bitset<5*maxK> kmerAminoXX::mask = mask.set();

/**
 * This function initializes the k-mer length and bit-mask.
 *
 * @param kmer_length k-mer length
 * @param bins the number of hash tables
 */
void kmerAminoXX::init(uint64_t& kmer_length, uint64_t& bins) {
    k = kmer_length; mask.reset();
    for (uint64_t i = 0; i < 5*k; ++i) {
        mask <<= 01u;    // fill all bits within the k-mer length with ones
        mask |= 01u;    // the remaining zero bits can be used to mask bits
    }

    bin = 0; // init the carry for the empty kmer
    table_count = bins; // the number of hash tables to use
    uint64_t last = 1 % table_count;
    for (int i = 1; i <= 2*(k + 1); i++)
    {
	    period.push_back(last);
	    last = (2 * last) % bins;
    }
}

/**
 * This function shifts a k-mer adding a new character to the left.
 * * http://www.cplusplus.com/forum/general/97378/
 * @param kmer bit sequence
 * @param c left character
 * @return right character
 */
char kmerAminoXX::shift_left(bitset<5 * maxK>& kmer, char& c) {
    uint64_t left = util::amino_char_to_bits(c);    // new leftmost character
    uint64_t right = 16*kmer[4]+8*kmer[3]+4*kmer[2]+2*kmer[1]+1*kmer[0];    // old rightmost character

    kmer >>= 05u;    // shift all current bits to the right by five positions

    for(int i = 4; i>=0; i--){   // encode the new character within the leftmost five bits
        kmer[5*k-(5-i)] = (left >> i) & 0x1;
    }

    kmer &= mask;    // set all bits to zero that exceed the k-mer length

    return util::amino_bits_to_char(right);    // return the dropped rightmost character
}


/**
 * This function shifts a k-mer adding a new character to the right.
 * http://www.cplusplus.com/forum/general/97378/
 * @param kmer bit sequence
 * @param c right character
 * @return left character
 */
char kmerAminoXX::shift_right(bitset<5 * maxK>& kmer, char& c) {
    uint64_t left = 16*kmer[5*k-1]+8*kmer[5*k-2]+4*kmer[5*k-3]+2*kmer[5*k-4]+1*kmer[5*k-5];    // old leftmost character
    uint64_t right = util::amino_char_to_bits(c);    // new rightmost character

    kmer <<= 05u;    // shift all current bits to the left by five positions
    for(int i = 4; i>=0; i--){// encode the new character within the rightmost five bits
        kmer[i] = (right >> i) & 0x1;
    }

    kmer &= mask;    // set all bits to zero that exceed the k-mer length

    // update the binning carry (solution of the shift-update-carry equation)
    bin = 32 * table_count // bias
		+ 16 * bin // Shift
		- 16 * period[2*k-1] * kmer[2*k-1] - 8 * kmer[2*k-2] * period[2*k - 2] - 4 * kmer[2*k-3] * period[2*k - 3]
        - 2 * kmer[2*k-4] * period[2*k - 4] - kmer[2*k-4] * period[2*k - 4];

    // Update
    for(int i = 0; i<=4; i++){
        bin += period[i] & ((right >> i) & 0x1);
    }

    bin %= table_count;

    return util::amino_bits_to_char(left);    // return the dropped leftmost character
}

/**
 * This function constructs the reverse complement of a given k-mer.
 *
 * @param kmer bit sequence
 * @param minimize only invert, if smaller
 * @return 1 if inverted, 0 otherwise
 */
bool kmerAminoXX::reverse_complement(bitset<5 * maxK>& kmer, bool minimize) {
    return false;
}


/**
 * The function computes the module on the Amino bitset
 *
 */
uint64_t kmerAminoXX::bit_mod(const bitset<5 * maxK>& kmer, uint64_t& module) {
	bitset<5 * maxK> bits = kmer;
	if (module <= 1){return 0;}
	
	uint64_t carry = 1;
	uint64_t rest = 0;

	if (bits[0]){rest++;} // Test the last bit
	
	for (uint64_t it=1; it < 5*k; it--){
	    carry = (2*carry) % module;
	    if (bits[it]){rest += carry;}
	}
	return rest % module;
}

