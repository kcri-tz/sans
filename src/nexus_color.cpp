#include "nexus_color.h"

/**
 * Plan:
 * 1. Create nexus file with necessary data TODO Add confidences
 * 2. Run splitstree to add network info TODO splitstree problems if too many splits
 * 3. Modify latest file with colors (see B) )
 */

/**
 * Help
 * - Set up project (in clion, makefile)
 * - Makefile correct?
 * - How to find the path to splitstree / require in PATH?
 * - Splitstree: zu viele splits werden nicht gezeichnet
 *      Warnung ausgegeben
 *      Mgl für Filter einbauen, der nur für nexus datei angewandt wird?
 * - Export Graphics from splitstree? -> -pdf -p
*/

// nexus datei ohne network erstellen und das nur erstellen lassen, wenn gewünscht,
// oder wenn pdf gewünscht

//für pdf -p -pdf (dann auch nur pdf, nexus wieder löschen)

// 2 dateien: taxa -> kategorie, kategorie -> farbe
// beides nach einem parameter (-l, --label oder so)
// farbliche labels gehen nur wenn pdf oder nexus gewünscht
// evtl nexus nur das einfache, wenn nexus und pdf gewünscht


/*
 * TODO at the moment:
 *  B) translate python script to add color (pimpnexus)
 *  D) check if makefile is correct
 *  E) Splitstree: zu viele splits werden nicht gezeichnet
 *      Bei nexus Ausgabe automatisch (-f strict fordern) #splits limitieren?
 *      **Bsp an Roland**
 *
 *  A) !check if working!
 *  open my nexus output with splitstree
 *  save modified file
 *  C) !check if complete!
 *  only have nexus output if desired (add flag nexus = true)
 *
 * num = denom_names.size() != denom_file_count (können mehrere files für eine sequenz sein)
 * DOCH s. main 683
 *
 * pimpnexus.py
 * ~/sans/example_data$ python3 ../scripts/pimpnexus.py ./test.nex ./color_file  > ./testpimped.nex
 *
 * sans2nexus.px
 * sans2conf_nexus.py
 *  -> add bootstrapping confidence values to nexus file CONFIDENCE=YES
 **/

bool program_in_path(const string& programName) {
    string command = "which " + programName;
    // Run command and check result
    int result = system(command.c_str());
    // If result is 0, the program is in PATH
    return result == 0;
}

void nexus_color::mod_via_splitstree(const string& nexus_file, const string& pdf, bool verbose, const string splitstree_path){
    // = "../splitstree4/SplitsTree"
    // TODO give option to give path tp splitstree

    // temporary file for splitstree commands
    const char* temp = "./temp_splitstree_commands";
    ofstream temp_file(temp);
    if (temp_file.is_open()) {

        // Writing commands for splitstree TODO full path
        temp_file << "begin SplitsTree;\nEXECUTE FILE=" << nexus_file << endl;
        temp_file << "UPDATE\nSAVE FILE=" << nexus_file << " REPLACE=YES\n";
        if(!pdf.empty()){ // TODO check if filename existent and replace or not
            temp_file << "EXPORTGRAPHICS format=PDF file=" << pdf << " REPLACE=yes\n";
        }
        temp_file << "QUIT\nend;";
        temp_file.close();

        // Running splitstree
        if(verbose) cout  << "Running SplitsTree\n";
        /// .../SplitsTree -g -S -c run_splitstree
        string command = splitstree_path + " -g -S -c " + temp;
        const char* splitstree_command = command.c_str();
        int result = system(splitstree_command); // executing command
        if(result == 0){
            if(verbose) cout << "Network created\n";
        } else {
            cerr << "Error while running Splitstree " << splitstree_path << endl;
            if (!program_in_path(splitstree_path)) {
                cerr << splitstree_path << " is not in the PATH." << endl;
            }
        }
        remove(temp);
    } else { // TODO
        cerr << "Unable to create temp command file for splitstree: " << temp << "\nTherefor unable to add network to nexus file" << endl;
    }
}


/// in construction ///
void nexus_color::color_nexus(string color_file, string nexus_file){
    // TODO ask for color file/families in the beginning? (so adding another variable?)
    //  separate script?
    bool translate = false;
    bool vertices = false;

    ifstream plain_nexus(nexus_file);
    ofstream colored_nexus("col_"+nexus_file);

    if(!plain_nexus.is_open()){
        cerr << "Error opening nexus file to modify with color.\n";
    }
    if(!colored_nexus.is_open()){
        cerr << "Error creating colored nexus file.\n";
    }

    //vector<string> lines;
    string line;
    // TODO read and save color mapping
    // read nexus data
    while (getline(plain_nexus, line)) {

        // BEGIN network; ... TRANSLATE [number][taxname]
        // -> match color via taxname to number
        // -> VERTICES [number] [vertice] <bg=...> <fg=...>

        // TODO case-insensitive search
        // identify block
        if(line.find("TRANSLATE") != string::npos){// starting point of: [number] [taxname]
            translate = true;
        }
        if (line.find("VERTICES") != string::npos) {// starting point of vertices
            vertices = true;
        }
        if(line.find(";") != string::npos){// end of any block
            translate = false;
            vertices = false;
        }

        // work with found blocks
        if(translate){
            // TODO match color -> number
        }
        if(vertices){
            // adding color
            //string color = "bg="+bg1+" "+bg2+" "+bg3 +" fg="+fg1+" "+fg2+" "+fg3+",\n";
            //line += color;
        }

        colored_nexus << line;
    }

    plain_nexus.close();
    colored_nexus.close();


    // receive file taxa <-> color
    // read file nexus with already added network via splitstree
    // map color <-> taxa <-> node
    // add color at respective nodes

    // special treatment if nodes would have multiple colors
    //      (draw two dots? other shape? different inner/outer color?)
    // nice-to-have: color for colorblind, or diff shape
}