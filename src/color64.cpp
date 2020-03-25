#include "color64.h"

/**
 * This is the number of colors.
 */
uint64_t color64::n;

/**
 * This is a bit-mask to erase all bits that exceed the color number.
 */
uint64_t color64::mask;

/**
 * This function initializes the color number and bit-mask.
 *
 * @param color_number color number
 */
void color64::init(uint64_t& color_number) {

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
void color64::set(uint64_t& color, uint64_t& pos) {
    color |= (0b1ull << pos) & mask;
}

/**
 * This function sets the bit at the given position to false.
 *
 * @param color bit sequence
 * @param pos position
 */
void color64::erase(uint64_t& color, uint64_t& pos) {
    color &= ~(0b1ull << pos) & mask;
}

/**
 * This function tests if the bit at the given position is set.
 *
 * @param color bit sequence
 * @param pos position
 * @return 1 if bit is set, 0 otherwise
 */
bool color64::test(uint64_t& color, uint64_t& pos) {
    return (color >> pos) & 0b1u;
}

/**
 * This function inverts a color set to reduce the number of ones.
 *
 * @param color bit sequence
 * @param minimize only invert, if smaller
 * @return 1 if inverted, 0 otherwise
 */
bool color64::complement(uint64_t& color, bool minimize) {

    uint64_t bits = color;    // copy the original color
    uint64_t ones = 0;    // counter for the number of ones

    for (uint64_t i = 0; i < n; ++i) {
        ones += bits & 0b1u;    // count the last bit
        bits >>= 01u;    // shift to the next bit pos.
    }

    // if minimize == true, return the color set with fewer ones
    if (minimize && (2*ones < n || (2*ones == n && (color & 0b1u)))) {
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
bool color64::is_subset(uint64_t& c1, uint64_t& c2) {
    return (c1 & c2) == c1;
}

/**
 * This function tests if two colors have no elements in common.
 *
 * @param c1 bit sequence
 * @param c2 bit sequence
 * @return true, if c1 & c2 are disjoint
 */
bool color64::is_disjoint(uint64_t& c1, uint64_t& c2) {
    return (c1 & c2) == 0b0u;
}

/**
 * This function tests if two colors are either subsets or disjoint.
 *
 * @param c1 bit sequence
 * @param c2 bis sequence
 * @return true, if c1 & c2 are compatible
 */
bool color64::is_compatible(uint64_t& c1, uint64_t& c2) {
    uint64_t intersect = c1 & c2;
    return intersect == 0b0u || intersect == c1 || intersect == c2;
}
