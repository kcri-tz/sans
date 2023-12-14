#include "nexus_color.h"

/**
 * Plan:
 * 1. Create nexus file with necessary data TODO Add confidences
 * 2. Run splitstree to add network info
 * 3. Modify latest file with colors
 */

/*
 * TODO at the moment:
 *  A) !check if working!
 *  open my nexus output with splitstree
 *  save modified file
 *  B)
 *  translate python script to add color (pimpnexus)
 *
 *  C) !check if complete!
 *  only have nexus output if desired (add flag nexus = true)
 *  D)
 *  check if makefile is correct
 *
 *   * TODO
 * - What about the output if only a newick/nexus output is given?
 *      Output should not be in sans at all, right? ja
 * - num = denom_names.size() != denom_file_count (können mehrere files für eine sequenz sein)
 * - Splitstree: zu viele splits werden nicht gezeichnet
 *      Bei nexus Ausgabe automatisch (-f strict fordern) #splits limitieren?
 *      bsp an Roland
 */

/* pimpnexus.py
 * ~/sans/example_data$ python3 ../scripts/pimpnexus.py ./test.nex ./color_file  > ./testpimped.nex
 *
 * sans2nexus.px
 * sans2conf_nexus.py
 *  -> add bootstrapping confidence values to nexus file CONFIDENCE=YES
 *
 * main.cpp 1092
 **/

void nexus_color::mod_via_splitstree(const string& nexus_file, bool verbose, const string splitstree_path){
    // = "../splitstree4/SplitsTree"
    // TODO adjust paths of Splitstree, nexus_file !!
    //  if verbose: print some info

    // temporary file for splitstree commands
    const char* temp = "./temp_splitstree_commands";
    ofstream temp_file(temp);
    if (temp_file.is_open()) {

        // Writing commands for splitstree TODO full path
        temp_file << "begin SplitsTree;\nEXECUTE FILE=" << nexus_file << endl;
        temp_file << "UPDATE\nSAVE FILE=" << nexus_file << " REPLACE=YES\n";
        //temp_file << "EXPORTGRAPHICS format=JPG file=" << image_file << "REPLACE=yes\n";
        temp_file << "QUIT\nend;";
        temp_file.close();

        // Running splitstree
        if(verbose) cout  << "Running SplitsTree\n";
        /// .../SplitsTree -g -S -c run_splitstree
        string command = splitstree_path + " -g -S -c " + temp;
        const char* splitstree_command = command.c_str();
        int result = system(splitstree_command); // executing command
        if(result == 0){
            if(verbose) cout << "Added network to nexus file\n";
        } else {
            cerr << "Could not run Splitstree " << splitstree_path << endl;
        }
        remove(temp);
    } else { // TODO
        cerr << "Unable to create temp command file for splitstree: " << temp << "\nTherefor unable to add network to nexus file" << endl;
    }
}

void nexus_color::color_nexus(string color_file, string nexus_file){
    // receive file taxa <-> color
    // read file nexus with already added network via splitstree
    // map color <-> taxa <-> node
    // add color at respective nodes
    // special treatment if nodes would have multiple colors
    //      (draw two dots? other shape? different inner/outer color?)
    // nice-to-have: color for colorblind, or diff shape
}

/**
 * - Coloring doesn't work properly ?
 * Help
 * - Set up project in clion
 * - Including nexus_color properly in the main, adding to makefile (?)
 * - Makefile correct?
*/