#include "Splits.h"

struct Trie {

    unordered_map<string, Trie *> nodes; // list of children
    unsigned int value = 0; // weight of the node
    unsigned int inverse = 0; // inverse weight
};

void addColors(Trie *trie, ColoredCDBG<> &graph,
               UnitigColors::const_iterator it, UnitigColors::const_iterator end, unsigned int weight) {

    vector<string> graph_colors = graph.getColorNames();
    vector<string> split_colors; // lists of colors in the graph and the split

    auto num = distance(it, end); // number of colors at the unitig
    auto max = graph_colors.size(); // number of all colors in the graph
    bool is_not_inverse;

    if (num == 0 || num == max) { // ignore empty or full splits
        return;
    } else if (2 * num == max) { // check the size of the split
        if (graph_colors.front() == graph.getColorName(it.getColorID())) {
            ++max;
        }
    }

    if (2 * num < max) { // collect the colors of the split
        is_not_inverse = true;
        for (; it != end; ++it) {
            split_colors.push_back(graph.getColorName(it.getColorID()));
        }
    } else { // collect the colors of the inverse split
        is_not_inverse = false;
        string current = graph.getColorName(it.getColorID());
        for (string &color : graph_colors) {
            if (it == end || color != current) {
                split_colors.push_back(color);
            } else { // filter out the colors at the unitig
                current = graph.getColorName((++it).getColorID());
            }
        }
    }

    Trie *subtrie = trie; // create a new subtrie for the split
    for (string &color : split_colors) {
        if (subtrie->nodes.count(color) == 0) {
            subtrie->nodes[color] = new Trie();
        } // follow the path down the trie and add the colors
        subtrie = subtrie->nodes[color];
    }
    if (is_not_inverse) { // store the weight of the split
        subtrie->value += weight;
    } else { // store the weight of the inverse split
        subtrie->inverse += weight;
    }
}

void getSplits(Trie *trie, string path, multimap<double, string, greater<double>> *splits) {

    if (trie->value > 0 || trie->inverse > 0) { // check if the split was observed
        //auto weight = trie->value + trie->inverse;
        auto weight = sqrt(trie->value) * sqrt(trie->inverse);
        splits->emplace(weight, path); // combine the weight and inverse weight
    }

    while (!trie->nodes.empty()) { // search colors depth first in the child node
        auto node = trie->nodes.begin();
        getSplits(node->second, path + "\t" + node->first, splits);
        trie->nodes.erase(node); // destroy the nodes that were visited
    }

    delete trie; // free the memory
}

void printFile(Trie *trie, ostream &out) {

    auto *splits = new multimap<double, string, greater<double>>();
    getSplits(trie, "", splits); // store the splits in a list sorted by weight

    while (!splits->empty()) { // print the weights and colors out to a file
        auto split = splits->begin();
        out << split->first << split->second << endl;
        splits->erase(split); // destroy the splits that were visited
    }

    delete splits; // free the memory
}

void buildTrie(ColoredCDBG<> &graph, string file_name) {

    Trie *trie = new Trie(); // create a new empty trie
    ofstream file_out(file_name); // open the file for output
    ostream out(file_out.rdbuf());

    for (const auto &unitig : graph) { // for each unitig of the graph

        const size_t nb_km = unitig.size - Kmer::k + 1; // number of k-mers in the unitig
        UnitigColors *uc_kmers = new UnitigColors[nb_km]; // create one color set per k-mer

        // get all pairs (k-mer position, color) for current unitig
        const UnitigColors *uc_unitig = unitig.getData()->getUnitigColors(unitig);
        const UnitigMapBase um(0, 1, Kmer::k, true);

        // iterator for pairs (k-mer position, color) for current unitig
        UnitigColors::const_iterator it = uc_unitig->begin(unitig);
        UnitigColors::const_iterator it_end = uc_unitig->end();

        for (; it != it_end; ++it) { // update k-mer color sets from iterators
            uc_kmers[it.getKmerPosition()].add(um, it.getColorID());
        }

        size_t len_segment = 1;
        //size_t len_segment = Kmer::k;

        // for each k-mer position in the current unitig
        for (size_t i = 1; i != nb_km; ++i, ++len_segment) {

            // there is a split (difference of color sets) if the color sets have different sizes
            bool split = (uc_kmers[i].size(um) != uc_kmers[i - 1].size(um));

            if (!split) {

                UnitigColors::const_iterator curr_it = uc_kmers[i].begin(um);
                UnitigColors::const_iterator curr_it_end = uc_kmers[i].end();

                UnitigColors::const_iterator prev_it = uc_kmers[i - 1].begin(um);
                UnitigColors::const_iterator prev_it_end = uc_kmers[i - 1].end();

                for (; !split && (curr_it != curr_it_end) &&
                       (prev_it != prev_it_end); ++curr_it, ++prev_it) {

                    // there is a split (difference of color sets) if the colors do not match
                    // between two overlapping k-mers
                    if (curr_it.getColorID() != prev_it.getColorID()) split = true;
                }
            }

            if (split) { // if there is a split

                addColors(trie, graph, uc_kmers[i - 1].begin(um), uc_kmers[i - 1].end(), len_segment);
                len_segment = 0; // add the colors to the trie and store the segment length
            }
        }

        addColors(trie, graph, uc_kmers[nb_km - 1].begin(um), uc_kmers[nb_km - 1].end(), len_segment);
        len_segment = 0; // add the colors to the trie and store the segment length

        delete[] uc_kmers; // free the memory
    }

    printFile(trie, out); // output the weights and colors of the splits to a file
}
