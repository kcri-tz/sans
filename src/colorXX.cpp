#include "colorXX.h"

/**
 * This is the number of colors.
 */
uint64_t colorXX::n;

/**
 * This is a bit-mask to erase all bits that exceed the color number.
 */
bitset<maxN> colorXX::mask;

/**
 * This function initializes the color number and bit-mask.
 *
 * @param color_number color number
 */
void colorXX::init(uint64_t& color_number) {

    n = color_number;
    for (uint64_t i = 0; i < n; ++i) {
        mask <<= 01u;    // fill all bits within the color number with ones
        mask |= 01u;    // the remaining zero bits can be used to mask bits
    }
}

/**
 * This function sets the bit at the given position to true.
 *
 * @param color bit sequence
 * @param pos position
 */
void colorXX::set(bitset<maxN>& color, uint64_t& pos) {
    color.set(pos) &= mask;
}

/**
 * This function sets the bit at the given position to false.
 *
 * @param color bit sequence
 * @param pos position
 */
void colorXX::erase(bitset<maxN>& color, uint64_t& pos) {
    color.reset(pos) &= mask;
}

/**
 * This function tests if the bit at the given position is set.
 *
 * @param color bit sequence
 * @param pos position
 * @return 1 if bit is set, 0 otherwise
 */
bool colorXX::test(bitset<maxN>& color, uint64_t& pos) {
    return color.test(pos);
}

/**
 * This function inverts a color set to reduce the number of ones.
 *
 * @param color bit sequence
 * @param minimize only invert, if smaller
 * @return 1 if inverted, 0 otherwise
 */
bool colorXX::complement(bitset<maxN>& color, bool minimize) {

    uint64_t ones = color.count();    // count the number of ones

    // if minimize == true, return the color set with fewer ones
    if (minimize && (2*ones < n || (2*ones == n && color[0]))) {
        return 0;    // not inverted
    } else {
        color = ~color & mask;    // flip the bits
        return 1;    // inverted
    }
}

/**
 * This function tests if one color set is a subset of another.
 *
 * @param c1 bit sequence
 * @param c2 bit sequence
 * @return true, if c1 is a subset of c2
 */
bool colorXX::is_subset(bitset<maxN>& c1, bitset<maxN>& c2) {
    return (c1 & c2) == c1;
}

/**
 * This function tests if two colors have no elements in common.
 *
 * @param c1 bit sequence
 * @param c2 bit sequence
 * @return true, if c1 & c2 are disjoint
 */
bool colorXX::is_disjoint(bitset<maxN>& c1, bitset<maxN>& c2) {
    return (c1 & c2) == 0b0u;
}

/**
 * This function tests if two colors are either subsets or disjoint.
 *
 * @param c1 bit sequence
 * @param c2 bis sequence
 * @return true, if c1 & c2 are compatible
 */
bool colorXX::is_compatible(bitset<maxN>& c1, bitset<maxN>& c2) {
    bitset<maxN> intersect = c1 & c2;
    return intersect == 0b0u || intersect == c1 || intersect == c2;
}
