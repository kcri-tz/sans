#include <iostream>

using namespace std;

/**
 * This class contains functions for working with color sets (maxN <= 64).
 */
class color64 {

private:

    /**
     * This is a bit-mask to erase all bits that exceed the color number.
     */
    static uint64_t mask;

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
    static void set(uint64_t& color, uint64_t& pos);

    /**
     * This function sets the bit at the given position to false.
     *
     * @param color bit sequence
     * @param pos position
     */
    static void erase(uint64_t& color, uint64_t& pos);

    /**
     * This function tests if the bit at the given position is set.
     *
     * @param color bit sequence
     * @param pos position
     * @return 1 if bit is set, 0 otherwise
     */
    static bool test(uint64_t& color, uint64_t& pos);

    /**
     * This function inverts a color set to reduce the number of ones.
     *
     * @param color bit sequence
     * @param minimize only invert, if smaller
     * @return 1 if inverted, 0 otherwise
     */
    static bool complement(uint64_t& color, bool minimize);

    /**
     * This function tests if one color set is a subset of another.
     *
     * @param c1 bit sequence
     * @param c2 bit sequence
     * @return true, if c1 is a subset of c2
     */
    static bool is_subset(uint64_t& c1, uint64_t& c2);

    /**
     * This function tests if two colors have no elements in common.
     *
     * @param c1 bit sequence
     * @param c2 bit sequence
     * @return true, if c1 & c2 are disjoint
     */
    static bool is_disjoint(uint64_t& c1, uint64_t& c2);

    /**
     * This function tests if two colors are either subsets or disjoint.
     *
     * @param c1 bit sequence
     * @param c2 bis sequence
     * @return true, if c1 & c2 are compatible
     */
    static bool is_compatible(uint64_t& c1, uint64_t& c2);

protected:

};
