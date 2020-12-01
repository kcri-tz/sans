#include "graphAmino.h"

/**
 * This is the size of the top list.
 */
uint64_t graphAmino::t;

/**
 * This is a hash table mapping k-mers to colors [O(1)].
 */
hash_map<kmerAmino_t, color_t> graphAmino::kmer_table;

/**
 * This is a hash table mapping colors to weights [O(1)].
 */
hash_map<color_t, array<uint32_t,2>> graphAmino::color_table;

/**
 * This is an ordered tree collecting the splits [O(log n)].
 */
multimap<double, color_t, greater<>> graphAmino::split_list;

/**
* These are the allowed chars.
*/
vector<char> graphAmino::allowedChars;

/**
 * This is a comparison function extending std::bitset.
 */
#if maxK > 12 || maxN > 64
namespace std {
    template <uint64_t N>
    bool operator<(const bitset<N>& x, const bitset<N>& y) {
        for (uint64_t i = N-1; i != -1; --i) {
            if (x[i] ^ y[i]) return y[i];
        }
        return false;
    }
}
#endif

/**
 * Initializes a new node struct.
 *
 * @param taxa color_t coding all taxa beneath this node
 * @param subsets list of subsets
 */
struct aminonode* newSet(color_t taxa, double weight, vector<aminonode*> subsets) {
    // declare and allocate new aminonode
    auto* aminonode = new struct aminonode();
    aminonode->taxa = std::move(taxa);
    aminonode->weight = std::move(weight);
    aminonode->subsets = std::move(subsets);
    return(aminonode);
}

/**
 * This function initializes the top list size and the allowed chars.
 *
 * @param t top list size
 */
void graphAmino::init(uint64_t& top_size) {
    t = top_size;

    graphAmino::allowedChars.push_back('A');
    graphAmino::allowedChars.push_back('B');
    graphAmino::allowedChars.push_back('C');
    graphAmino::allowedChars.push_back('D');
    graphAmino::allowedChars.push_back('E');
    graphAmino::allowedChars.push_back('F');
    graphAmino::allowedChars.push_back('G');
    graphAmino::allowedChars.push_back('H');
    graphAmino::allowedChars.push_back('I');
    graphAmino::allowedChars.push_back('J');
    graphAmino::allowedChars.push_back('K');
    graphAmino::allowedChars.push_back('L');
    graphAmino::allowedChars.push_back('M');
    graphAmino::allowedChars.push_back('N');
    graphAmino::allowedChars.push_back('O');
    graphAmino::allowedChars.push_back('P');
    graphAmino::allowedChars.push_back('Q');
    graphAmino::allowedChars.push_back('R');
    graphAmino::allowedChars.push_back('S');
    graphAmino::allowedChars.push_back('T');
    graphAmino::allowedChars.push_back('U');
    graphAmino::allowedChars.push_back('V');
    graphAmino::allowedChars.push_back('W');
    graphAmino::allowedChars.push_back('X');
    graphAmino::allowedChars.push_back('Y');
    graphAmino::allowedChars.push_back('Z');
    graphAmino::allowedChars.push_back('*');

}

/**
 * This function extracts k-mers from a sequence and adds them to the hash table.
 *
 * @param str amino acid sequence
 * @param color color flag
 * @param reverse merge complements
 */
void graphAmino::add_kmers(string& str, uint64_t& color, bool& reverse) {
    if (str.length() < kmerAmino::k) return;    // not enough characters

    uint64_t pos;    // current position in the string, from 0 to length
    kmerAmino_t kmerAmino;    // create a new empty bit sequence for the k-mer
    kmerAmino_t rcmer;    // create a bit sequence for the reverse complement

    uint64_t begin = 0;
next_kmer:
    pos = begin;

    for (; pos < str.length(); ++pos) {    // collect the bases from the string
        if (!isAllowedChar(pos, str)) {
            begin = pos+1;    // str = str.substr(pos+1, string::npos);
            goto next_kmer;    // unknown base, start a new k-mer from the beginning
        }
        kmerAmino::shift_right(kmerAmino, str[pos]);    // shift each base into the bit sequence

        if (pos+1 - begin >= kmerAmino::k) {
            rcmer = kmerAmino;
            if (reverse) kmerAmino::reverse_complement(rcmer, true);    // invert the k-mer, if necessary
            color::set(kmer_table[rcmer], color);    // update the k-mer with the current color
        }
    }
}

/**
 * This function extracts k-mer minimizers from a sequence and adds them to the hash table.
 *
 * @param str amino acid sequence
 * @param color color flag
 * @param reverse merge complements
 * @param m number of k-mers to minimize
 */
void graphAmino::add_minimizers(string& str, uint64_t& color, bool& reverse, uint64_t& m) {
    if (str.length() < kmerAmino::k) return;    // not enough characters

    vector<kmerAmino_t> sequence_order;    // k-mers ordered by their position in sequence
    multiset<kmerAmino_t> value_order;    // k-mers ordered by their lexicographical value

    uint64_t pos;    // current position in the string, from 0 to length
    kmerAmino_t kmerAmino;    // create a new empty bit sequence for the k-mer
    kmerAmino_t rcmer;    // create a bit sequence for the reverse complement

    uint64_t begin = 0;
next_kmer:
    pos = begin;
    sequence_order.clear();
    value_order.clear();

    for (; pos < str.length(); ++pos) {    // collect the bases from the string
        if (!isAllowedChar(pos, str)) {
            begin = pos+1;    // str = str.substr(pos+1, string::npos);
            goto next_kmer;    // unknown base, start a new k-mer from the beginning
        }
        kmerAmino::shift_right(kmerAmino, str[pos]);    // shift each base into the bit sequence

        if (pos+1 - begin >= kmerAmino::k) {
            rcmer = kmerAmino;
            if (reverse) kmerAmino::reverse_complement(rcmer, true);    // invert the k-mer, if necessary

            if (sequence_order.size() == m) {
                value_order.erase(*sequence_order.begin());    // remove k-mer outside the window
                sequence_order.erase(sequence_order.begin());
            }
            value_order.emplace(rcmer);    // insert k-mer ordered by its lexicographical value
            sequence_order.emplace_back(rcmer);

            if (sequence_order.size() == m) {
                color::set(kmer_table[*value_order.begin()], color);    // update the k-mer with the current color
            }
        }
    }
}
/**
 * This function checks if the character at the given position is allowed.
 * @param pos position in str
 * @param str the current part of the sequence
 * @return true if allowed, false otherwise
 */
bool graphAmino::isAllowedChar(uint64_t pos, string &str) {
    bool allowed = false;
    char &currentChar = str[pos];

    for(int i = 0; i<graphAmino::allowedChars.size() && !allowed; i++){
        allowed =  graphAmino::allowedChars.at(i) == currentChar;
    }

    return allowed;
}

/**
 * This function iterates over the hash table and calculates the split weights.
 *
 * @param mean weight function
 * @param verbose print progress
 */
void graphAmino::add_weights(double mean(uint32_t&, uint32_t&), bool& verbose) {
    double min_value = numeric_limits<double>::min();    // current min. weight in the top list (>0)
    uint64_t cur = 0, prog = 0, next;
    uint64_t max = kmer_table.size();
loop:
    for (auto it = kmer_table.begin(); it != kmer_table.end(); ++it) {    // iterate over k-mer hash table
        if (verbose) {
            next = 100*cur/max;
             if (prog < next)  cout << "\33[2K\r" << "Processing splits... " << next << "%" << flush;
            prog = next; cur++;
        }
        color_t& color = it.value();    // get the color set for each k-mer
        bool pos = color::complement(color, true);    // invert the color set, if necessary
        if (color == 0) continue;    // ignore empty splits
        array<uint32_t,2>& weight = color_table[color];    // get the weight and inverse weight for the color set

        double old_value = mean(weight[0], weight[1]);    // calculate the old mean value
        if (old_value >= min_value) {    // if it is greater than the min. value, find it in the top list
            auto range = split_list.equal_range(old_value);    // get all color sets with the given weight
            for (auto it = range.first; it != range.second; ++it) {
                if (it->second == color) {    // iterate over the color sets to find the correct one
                    split_list.erase(it);    // erase the entry with the old weight
                    break;
                }
            }
        }
        weight[pos]++;    // update the weight or the inverse weight of the current color set

        double new_value = mean(weight[0], weight[1]);    // calculate the new mean value
        if (new_value >= min_value) {    // if it is greater than the min. value, add it to the top list
            split_list.emplace(new_value, color);    // insert it at the correct position ordered by weight
            if (split_list.size() > t) {
                split_list.erase(--split_list.end());    // if the top list exceeds its limit, erase the last entry
                min_value = split_list.rbegin()->first;    // update the min. value for the next iteration
            }
        }
    }
}

/**
 * This function adds a single split (weight and colors) to the output list.
 *
 * @param weight split weight
 * @param color split colors
 */
void graphAmino::add_split(double& weight, color_t& color) {
    split_list.emplace(weight, color);    // insert it at the correct position ordered by weight
    if (split_list.size() > t) {
        split_list.erase(--split_list.end());    // if the top list exceeds its limit, erase the last entry
    }
}

/**
 * This function filters a greedy maximum weight tree compatible subset.
 *
 * @param verbose print progress
 */
void graphAmino::filter_strict(bool& verbose) {
    filter_strict(nullptr, verbose);
}

/**
 * This function filters a greedy maximum weight tree compatible subset and returns a newick string.
 *
 * @param map function that maps an integer to the original id, or null if no newick output wanted
 * @param verbose print progress
 */
string graphAmino::filter_strict(std::function<string(const uint64_t&)> map, bool& verbose) {
    auto tree = vector<color_t>();    // create a set for compatible splits
    auto it = split_list.begin();
    uint64_t cur = 0, prog = 0, next;
    uint64_t max = split_list.size();
loop:
    while (it != split_list.end()) {
        if (verbose) {
            next = 100*cur/max;
             if (prog < next)  cout << "\33[2K\r" << "Filtering splits... " << next << "%" << flush;
            prog = next; cur++;
        }
        if (test_strict(it->second, tree)) {
            tree.emplace_back(it->second);
            ++it; goto loop;    // if compatible, add the new split to the set
        }
        it = split_list.erase(it);    // otherwise, remove split
    }
    if (map) {
        aminonode* root = build_tree(tree);
        return print_tree(root, map) + ";\n";
    } else {
        return "";
    }
}

/**
 * This function filters a greedy maximum weight weakly compatible subset.
 *
 * @param verbose print progress
 */
void graphAmino::filter_weakly(bool& verbose) {
    auto network = vector<color_t>();    // create a set for compatible splits
    auto it = split_list.begin();
    uint64_t cur = 0, prog = 0, next;
    uint64_t max = split_list.size();
loop:
    while (it != split_list.end()) {
        if (verbose) {
            next = 100*(cur*cur)/(max*max);
             if (prog < next)  cout << "\33[2K\r" << "Filtering splits... " << next << "%" << flush;
            prog = next; cur++;
        }
        if (test_weakly(it->second, network)) {
            network.emplace_back(it->second);
            ++it; goto loop;    // if compatible, add the new split to the set
        }
        it = split_list.erase(it);    // otherwise, remove split
    }
}

/**
 * This function filters a greedy maximum weight n-tree compatible subset.
 *
 * @param n number of trees
 * @param verbose print progress
 */
void graphAmino::filter_n_tree(uint64_t n, bool& verbose) {
    filter_n_tree(n, nullptr, verbose);
}

/**
 * This function filters a greedy maximum weight n-tree compatible subset and returns a string with all trees in newick format.
 *
 * @param n number of trees
 * @param map function that maps an integer to the original id, or null
 * @param verbose print progress
 */
string graphAmino::filter_n_tree(uint64_t n, std::function<string(const uint64_t&)> map, bool& verbose) {
    auto forest = vector<vector<color_t>>(n);    // create a set for compatible splits
    auto it = split_list.begin();
    uint64_t cur = 0, prog = 0, next;
    uint64_t max = split_list.size();
loop:
    while (it != split_list.end()) {
        if (verbose) {
            next = 100*cur/max;
             if (prog < next)  cout << "\33[2K\r" << "Filtering splits... " << next << "%" << flush;
            prog = next; cur++;
        }
       for (auto& tree : forest)
        if (test_strict(it->second, tree)) {
            tree.emplace_back(it->second);
            ++it; goto loop;    // if compatible, add the new split to the set
        }
        it = split_list.erase(it);    // otherwise, remove split
    }
    // output
    string s;
    if (map) {
        for (auto& tree : forest) {
            aminonode* root = build_tree(tree);
            s += print_tree(root, map) + ";\n";
        }
    }
    return s;
}

/**
 * This function tests if a split is compatible with an existing set of splits.
 *
 * @param color new split
 * @param color_set set of splits
 * @return true, if compatible
 */
bool graphAmino::test_strict(color_t& color, vector<color_t>& color_set) {
    for (auto& elem : color_set) {
        if (!color::is_compatible(elem, color)) {
            return false;    // compare to each split in the set
        }
    }
    return true;
}

/**
 * This function tests if a split is weakly compatible with an existing set of splits.
 *
 * @param color new split
 * @param color_set set of splits
 * @return true, if weakly compatible
 */
bool graphAmino::test_weakly(color_t& color, vector<color_t>& color_set) {
    for (auto& elem1 : color_set) {
        for (auto& elem2 : color_set) {
            if (elem1 != elem2) {
                if (!color::is_weakly_compatible(elem1, elem2, color)) {
                    return false;    // compare to each split in the set
                }
            }
        }
    }
    return true;
}

/**
 * This function recursively refines a given set/tree structure by a given split.
 *
 * @param current_set node of currently considered (sub-)set/tree structure
 * @param split color set to refine by
 * @return whether or not the given split is compatible with the set/tree structure
 */
bool graphAmino::refine_tree(aminonode* current_set, color_t& split, color_t& allTaxa) {
    // possible cases:
    // splitsize <2: nothing has to be done
    // split equals one subset -> warning: split twice
    // split is fully contained in one subset -> recurse
    // inverse split ... (i.e. split covers one subset partially) -> recurse with inverse
    // split covers several subsets completely -> introduce new split
    if (color::size(split, false) < 2 || color::size(allTaxa, false) - color::size(split, false) < 2) { return true; }

    vector<aminonode*> *subsets = &current_set->subsets;
    vector<aminonode*> fullycoveredsubsets = {};
    aminonode* partiallycoveredsubset = nullptr;

    for (aminonode* subset : *subsets) {
        color_t subtaxa = subset->taxa;
        if (split == subtaxa) {
            return true;
        }
        // split.issubset(subtaxa)?
        if ((split & subtaxa) == split) { return refine_tree(subset, split, allTaxa); }
        // subtaxa.issubset(split):
        if ((subtaxa & split) == subtaxa) { fullycoveredsubsets.push_back(subset); }
        // elif not subtaxa.isdisjoint(split): # does intersect
        else if ((subtaxa & split) != 0b0u) {
            // if partiallycoveredsubset:
            if (partiallycoveredsubset != nullptr) { return false; } //there cannot be more than one
            else { partiallycoveredsubset = subset; }
        }
    }

    if (partiallycoveredsubset != nullptr) {
        if (fullycoveredsubsets.size() == subsets->size()-1) {
            // recurse into this subset with inverse split
            color_t inversesplit = color::complement(split, false);
            // if inversesplit.issubset(partiallycoveredsubset[1]):
            if ((inversesplit & partiallycoveredsubset->taxa) == inversesplit) {
                return refine_tree(partiallycoveredsubset, inversesplit, allTaxa);
            } else { return false; }
        } else { return false; }
    } else if (fullycoveredsubsets.size() > 1) {
        // introduce new split
        color_t newsubtaxa = 0b0u;
        for(aminonode* subset : fullycoveredsubsets) { newsubtaxa |= subset->taxa; }
        // get weight of split
        double weight = 0;
        auto it = split_list.begin();
        while (it != split_list.end()) {
            if (it->second == split){
                weight = it->first;
                break;
            }
            it++;
        }
        aminonode* newset = newSet(newsubtaxa, weight, fullycoveredsubsets);
        // remove old sets
        for(aminonode* subset : fullycoveredsubsets) {
            // subsets.remove(subset)
            subsets->erase(std::remove(subsets->begin(), subsets->end(), subset), subsets->end());
        }
        // add new set
        subsets->push_back(newset);
        return true;
    } else {
        std::cerr << "ERROR: this cannot be: just one fully covered subset and nothing else!?" << endl;
        exit(EXIT_FAILURE);
    }
}

/**
 * This function returns a tree structure (struct node) generated from the given list of color sets.
 *
 * @param color_set list of color sets
 * @return tree structure (struct node)
 */
aminonode* graphAmino::build_tree(vector<color_t>& color_set) {
    //initialize set of trivial splits
    vector<aminonode*> subsets = {};
    color_t allTaxa = 0b0u;

    for (uint64_t i = 0; i < color::n; i++) {
        color_t leaf = 0b0u;
        color::set(leaf, i);
        color::set(allTaxa, i);
        vector<aminonode*> emptyset = {};
        // get weight
        double weight = 0;
        auto it = split_list.begin();
        while (it != split_list.end()) {
            if (it->second == leaf){
                weight = it->first;
                break;
            }
            it++;
        }
        aminonode* newset = newSet(leaf, weight, emptyset);
        subsets.push_back(newset);
    }
    aminonode* sets = newSet(allTaxa, 0, subsets);

    for (color_t split : color_set) {
        // split if possible
        if (!refine_tree(sets, split, allTaxa)) {
            std::cerr << "ERROR: splits are incompatible" << endl;
            exit(EXIT_FAILURE);
        }
    }
    return sets;
}

/**
 * This function returns a newick string generated from the given tree structure (set).
 *
 * @param root root of the tree/set structure
 * @return newick string
 */
string graphAmino::print_tree(aminonode* root, std::function<string(const uint64_t&)> map) {
    vector<aminonode*> subsets = root->subsets;
    color_t taxa = root->taxa;

    if (subsets.empty()){    // leaf set
        if (color::size(taxa, false) == 0) {
            std::cerr << "ERROR: child with no taxon!?" << endl;
            exit(EXIT_FAILURE);
        } else if (color::size(taxa, false) == 1) {
            return map(color::pos(taxa)) + ":" + to_string(root->weight);
        } else {
            std::cerr << "ERROR: child with more than one taxon!?" << endl;
            exit(EXIT_FAILURE);
        }
    }
    else {
        string s = "(";
        for (aminonode* subset : subsets) {
            s += print_tree(subset, map);
            if (subset != subsets.back()) { s += ","; }
        }
        s += "):";
        s += to_string(root->weight);
        return s;
    }
}
