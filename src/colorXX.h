#include <iostream>
#include <bitset>

using namespace std;

/**
 * This is the max. number of colors defined as a preprocessor directive.
 */
#ifndef maxN
#define maxN 64
#endif

/**
 * This class contains functions for working with color sets (maxN > 64).
 */
class colorXX {

private:

    /**
     * This is a bit-mask to erase all bits that exceed the color number.
     */
    static bitset<maxN> mask;

public:

    /**
     * This is the number of colors.
     */
    static uint64_t n;

    /**
     * This function initializes the color number and bit-mask.
     *
     * @param color_number color number
     */
    static void init(uint64_t& color_number);

    /**
     * This function sets the bit at the given position to true.
     *
     * @param color bit sequence
     * @param pos position
     */
    static void set(bitset<maxN>& color, uint64_t& pos);

    /**
     * This function sets the bit at the given position to false.
     *
     * @param color bit sequence
     * @param pos position
     */
    static void erase(bitset<maxN>& color, uint64_t& pos);

    /**
     * This function tests if the bit at the given position is set.
     *
     * @param color bit sequence
     * @param pos position
     * @return 1 if bit is set, 0 otherwise
     */
    static bool test(bitset<maxN>& color, uint64_t& pos);

    /**
     * This function inverts a color set to reduce the number of ones.
     *
     * @param color bit sequence
     * @param minimize only invert, if smaller
     * @return 1 if inverted, 0 otherwise
     */
    static bool complement(bitset<maxN>& color, bool minimize);

    /**
     * This function tests if one color set is a subset of another.
     *
     * @param c1 bit sequence
     * @param c2 bit sequence
     * @return true, if c1 is a subset of c2
     */
    static bool is_subset(bitset<maxN>& c1, bitset<maxN>& c2);

    /**
     * This function tests if two colors have no elements in common.
     *
     * @param c1 bit sequence
     * @param c2 bit sequence
     * @return true, if c1 & c2 are disjoint
     */
    static bool is_disjoint(bitset<maxN>& c1, bitset<maxN>& c2);

    /**
     * This function tests if two colors are either subsets or disjoint.
     *
     * @param c1 bit sequence
     * @param c2 bis sequence
     * @return true, if c1 & c2 are compatible
     */
    static bool is_compatible(bitset<maxN>& c1, bitset<maxN>& c2);

protected:

};
