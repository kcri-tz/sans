#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cctype>
#include <unordered_map>

using namespace std;

#ifndef SRC_NEXUS_COLOR_H
#define SRC_NEXUS_COLOR_H

namespace nexus_color{

    /**
     * This function adds a network to the initially generated nexus file via SplitsTree.
     * @param nexus_file Path to the initial nexus file
     * @param verbose If info should be printed
     * @param splitstree_path Path to SplitsTree
     */
    void mod_via_splitstree(const string& nexus_file, const string& pdf, bool verbose = false, const string splitstree_path = "SplitsTree");

    void color_nexus(const string& nexus_file, const string& tax_fam_file, const string& fam_clr_file);
}

#endif //SRC_NEXUS_COLOR_H