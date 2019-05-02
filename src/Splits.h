#ifndef BIFROST_SPLITS_H
#define BIFROST_SPLITS_H

using namespace std;

/**
 * This data structure represents a trie of colors.
 * @param nodes list of children
 * @param value weight of the node
 * @param inverse inverse weight
 */
struct Trie;

/**
 * This function inserts the colors into the trie.
 * @param trie trie to construct
 * @param graph de bruijn graph
 * @param it unitig iterator begin
 * @param end unitig iterator end
 * @param weight weight of the node
 */
void addColors(Trie *trie, ColoredCDBG<> &graph,
               UnitigColors::const_iterator it, UnitigColors::const_iterator end, unsigned int weight);

/**
 * This function extracts the splits from the trie.
 * @param trie trie to traverse
 * @param path current trie path
 * @param splits list of splits
 */
void getSplits(Trie *trie, string path, multimap<double, string, greater<double>> *splits);

/**
 * This function prints the splits out to a file.
 * @param trie trie to traverse
 * @param out file output stream
 */
void printFile(Trie *trie, ostream &out);

/**
 * This function constructs a trie from a graph.
 * @param graph de bruijn graph
 * @param file_name output file
 */
void buildTrie(ColoredCDBG<> &graph, string file_name);

#endif //BIFROST_SPLITS_H
