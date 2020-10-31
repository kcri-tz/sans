#include "kmerAminoXX.h"

/**
 * This is the length of a k-mer.
 */
uint64_t kmerAminoXX::k;

/**
 * This is a bit-mask to erase all bits that exceed the k-mer length.
 */
bitset<5*maxK> kmerAminoXX::mask = mask.set();

/**
 * This function initializes the k-mer length and bit-mask.
 *
 * @param kmer_length k-mer length
 */
void kmerAminoXX::init(uint64_t& kmer_length) {
    k = kmer_length; mask.reset();
    for (uint64_t i = 0; i < 5*k; ++i) {
        mask <<= 01u;    // fill all bits within the k-mer length with ones
        mask |= 01u;    // the remaining zero bits can be used to mask bits
    }
}

/**
 * This function shifts a k-mer adding a new character to the left.
 * * http://www.cplusplus.com/forum/general/97378/
 * TODO debuggen!!!
 * @param kmer bit sequence
 * @param c left character
 * @return right character
 */
char kmerAminoXX::shift_left(bitset<5 * maxK>& kmer, char& c) {
    uint64_t left = char_to_bits(c);    // new leftmost character
    uint64_t right = 5*kmer[4]+kmer[3]+kmer[2]+kmer[1]+kmer[0];    // old rightmost character

    kmer >>= 05u;    // shift all current bits to the right by five positions

    for(int i = 4; i>=0; i--){   // encode the new character within the leftmost five bits
        kmer[5*k-(5-i)] = (left >> i) & 0x1;
    }

    kmer &= mask;    // set all bits to zero that exceed the k-mer length

    return bits_to_char(right);    // return the dropped rightmost character
}


/**
 * This function shifts a k-mer adding a new character to the right.
 * http://www.cplusplus.com/forum/general/97378/
 * TODO debuggen!!!
 * @param kmer bit sequence
 * @param c right character
 * @return left character
 */
char kmerAminoXX::shift_right(bitset<5 * maxK>& kmer, char& c) {
    uint64_t left = 5*kmer[5*k-1]+kmer[5*k-2]+kmer[5*k-3]+kmer[5*k-4]+kmer[5*k-5];    // old leftmost character
    uint64_t right = char_to_bits(c);    // new rightmost character

    kmer <<= 05u;    // shift all current bits to the left by five positions
    for(int i = 4; i>=0; i--){// encode the new character within the rightmost five bits
        kmer[i] = (right >> i) & 0x1;
    }

    kmer &= mask;    // set all bits to zero that exceed the k-mer length

    return bits_to_char(left);    // return the dropped leftmost character
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
 * This function encodes a single character to five bits.
 *
 * @param c character
 * @return bit sequence
 */
uint64_t kmerAminoXX::char_to_bits(char& c) {

    switch (c) {
        case 'P':
            return 0b00000u;
        case 'A':
            return 0b00001u;
        case 'G':
            return 0b00010u;
        case 'Q':
            return 0b00011u;
        case 'N':
            return 0b00100u;
        case 'E':
            return 0b00101u;
        case 'D':
            return 0b00110u;
        case 'T':
            return 0b00111u;
        case 'S':
            return 0b01000u;
        case 'C':
            return 0b01001u;
        case 'V':
            return 0b01010u;
        case 'I':
            return 0b01011u;
        case 'M':
            return 0b01100u;
        case 'L':
            return 0b01101u;
        case 'F':
            return 0b01110u;
        case 'Y':
            return 0b01111u;
        case 'W':
            return 0b10000u;
        case 'H':
            return 0b10001u;
        case 'K':
            return 0b10010u;
        case 'R':
            return 0b10011u;
        case 'X':
            return 0b10100u;
        case '*':
            return 0b10101u;
        default:
            cerr << "Error: Invalid character " << c << "." << endl;
            return -1;
    }
}

/**
 * This function decodes five bits to a single character.
 *
 * @param b bit sequence
 * @return character
 */
char kmerAminoXX::bits_to_char(uint64_t& b) {
    switch (b) {
        case 0b00000u:
            return 'P';
        case 0b00001u:
            return 'A';
        case 0b00010u:
            return 'G';
        case 0b00011u:
            return 'Q';
        case 0b00100u:
            return 'N';
        case 0b00101u:
            return 'E';
        case 0b00110u:
            return 'D';
        case 0b00111u:
            return 'T';
        case 0b01000u:
            return 'S';
        case 0b01001u:
            return 'C';
        case 0b01010u:
            return 'V';
        case 0b01011u:
            return 'I';
        case 0b01100u:
            return 'M';
        case 0b01101u:
            return 'L';
        case 0b01110u:
            return 'F';
        case 0b01111u:
            return 'Y';
        case 0b10000u:
            return 'W';
        case 0b10001u:
            return 'H';
        case 0b10010u:
            return 'K';
        case 0b10011u:
            return 'R';
        case 0b10100u:
            return 'X';
        case 0b10101u:
            return '*';
        default:
            cerr << "Error: Invalid bit sequence " << b << "." << endl;
            return -1;
    }
}
