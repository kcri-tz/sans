#include "nexus_color.h"

/**
 * Plan:
 * 1. Create nexus file with necessary data TODO Add confidences
 * 2. Run splitstree to add network info TODO splitstree problems if too many splits
 * 3. Modify latest file with colors (see B) )
 *
 * TODO at the moment:
 *  B) translate python script to add color (pimpnexus)
 *  A) add confidence values
 *
 * - Splitstree: zu viele splits werden nicht gezeichnet
 *    Warnung ausgegeben
 *    Mgl für Filter einbauen, der nur für nexus datei angewandt wird?
*/

// 2 dateien: taxa -> kategorie, kategorie -> farbe
// beides nach einem parameter (-l, --label oder so)
// farbliche labels gehen nur wenn pdf oder nexus gewünscht
// evtl nexus nur das einfache, wenn nexus und pdf gewünscht

/* num = denom_names.size() != denom_file_count (können mehrere files für eine sequenz sein)
 * DOCH s. main 683
 *
 * pimpnexus.py
 *
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

/** TODO !!TEST!!
 * Creates a map of a given tab separated file with two fields (key,value) per line.
 * Also creates map with the given values as keys for further usage.
 * @param map The dictionary consistent of the provided (key,value) pairs.
 * @param filename The input file to read from.
 * @param value_map The dictionary with the provided values as keys.
 */
void create_maps(unordered_map<string, string>& map, const string& filename, unordered_map<string, string>& value_map){

    string line;
    ifstream infile(filename);
    if(!infile.is_open()){
        cerr << "Error: Could not read input file " << filename << endl;
        return;
    }
    while (getline(infile, line)){
        // String stream to split and read line
        istringstream iss(line);
        string key, value; // two fields, tab separated
        if (getline(iss, key, '\t') && getline(iss, value)) {
            // Store fields in dictionary
            map[key] = value;
            value_map[value] = "";
        }
    }
    infile.close();

}

void nexus_color::mod_via_splitstree(const string& nexus_file, const string& pdf, bool verbose, const string splitstree_path){
    // = "../splitstree4/SplitsTree"
    // TODO give option to give path to splitstree

    // temporary file for splitstree commands
    const char* temp = "./temp_splitstree_commands";
    ofstream temp_file(temp);
    if (temp_file.is_open()) {

        // Writing commands for splitstree
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
void nexus_color::color_nexus(const string& nexus_file, const string& tax_fam_file, const string& fam_clr_file){
    // TODO ask for color file/families in the beginning
    // 2 dateien: taxa -> kategorie, kategorie -> farbe
    // beides nach einem parameter (-l, --label oder so)
    // TODO if only taxa -> family given, choose colors myself, otherwise use given colors
    bool translate = false;
    bool vertices = false;
    //vector<string> lines;
    string line;

    // TODO read and save tax/fam/color mapping
    //  (tax, clr) (via (tax, fam) -> (fam, clr))
    //  with color either given or self-calculated
    unordered_map<string, string> tax_fam_map;
    unordered_map<string, string> fam_clr_map;
    create_maps(tax_fam_map, tax_fam_file, fam_clr_map);
    // Reading and saving fam -> col mapping if given
    if(!fam_clr_file.empty()){
        // TODO add colors to fam_clr_map
        //  and check if provided fams stimmen überein
    } else {
        // TODO create colors according to #families
    }

    line = "";
    ifstream plain_nexus(nexus_file);
    ofstream colored_nexus("col_"+nexus_file);

    if(!plain_nexus.is_open()){
        cerr << "Error opening nexus file to modify with color.\n";
    }
    if(!colored_nexus.is_open()){
        cerr << "Error creating colored nexus file.\n";
    }

    // read nexus data
    while (getline(plain_nexus, line)) {

        // work with found blocks
        if(translate){
            // (number, taxname)
            // TODO match color -> number
            // BEGIN network; ... TRANSLATE [number][taxname]
            // -> match color via taxname to number
            // -> VERTICES [number] [vertice] <bg=...> <fg=...>

        }
        if(vertices){
            // adding color
            //string color = "bg="+bg1+" "+bg2+" "+bg3 +" fg="+fg1+" "+fg2+" "+fg3+",\n";
            //line += color;
        }

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