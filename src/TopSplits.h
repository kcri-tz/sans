#ifndef BIFROST_TOPSPLITS_H
#define BIFROST_TOPSPLITS_H

#include "TopTrie.cpp"

using namespace std;

void searchGraph(ColoredCDBG<> &graph, unsigned int top_num, string file_name);

void addColors(TopTrie *trie, ColoredCDBG<> &graph,
               UnitigColors::const_iterator it, UnitigColors::const_iterator end, unsigned int value);

void printSplits(TopTrie *trie, ostream &out);

#endif //BIFROST_TOPSPLITS_H
