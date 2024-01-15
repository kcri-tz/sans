#include "nexus_color.h"

struct rgb_color {
    int r, g, b;
};


bool nexus_color::program_in_path(const string& programName) {
    // suppress output of this command depending on system
    #ifdef _WIN32
    #define DEV_NULL "nul"
    #else
    #define DEV_NULL "/dev/null"
    #endif

    string command = "which " + programName + " > " + string(DEV_NULL);
    // Run command and check result
    int result = system(command.c_str());
    // If result is 0, the program is in PATH
    return result == 0;
}

/**
 * Creates a map of a given tab separated file with two fields (key,value) per line.
 * Also creates map with the given values as keys for further usage.
 * @param map The dictionary consistent of the provided (key,value) pairs.
 * @param filename The input file to read from.
 * @param value_map The dictionary with the provided values as keys.
 */
void create_maps(unordered_map<string, string>& map, const string& filename, unordered_map<string, rgb_color>& value_map){

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
            value_map[value]; // TODO add empty color?
        }
    }
    infile.close();
}

void nexus_color::mod_via_splitstree(const string& nexus_file, const string& pdf, bool verbose, const string splitstree_path){
    // = "../splitstree4/SplitsTree"
    if (!program_in_path(splitstree_path)) {
        cerr << splitstree_path << " is not in the PATH." << endl;
        return;
    }
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
        }
        remove(temp);
    } else { // TODO
        cerr << "Unable to create temp command file for splitstree: " << temp << "\nTherefor unable to add network to nexus file" << endl;
    }
}

/// in construction ///
void nexus_color::color_nexus(const string& nexus_file, const string& tax_grp_file, const string& grp_clr_file){

    bool translate = false;
    bool vertices = false;
    //vector<string> lines;
    string line;

    //  TODO self-calculated color

    unordered_map<string, string> tax_grp_map;
    unordered_map<string, rgb_color> grp_clr_map;
    // Reading and saving taxa -> grp mapping
    create_maps(tax_grp_map, tax_grp_file, grp_clr_map);

    // (Reading &) Saving grp -> col mapping
    if(!grp_clr_file.empty()){
        // add colors to grp_clr_map
        string grp_clr;
        ifstream infile(grp_clr_file);
        if(!infile.is_open()){
            cerr << "Error: Could not read input file " << grp_clr_file << endl;
            return;
        }
        while (getline(infile, grp_clr)){
            // String stream to split and read line
            istringstream iss(grp_clr);
            string grp, clr; // two fields, tab separated
            if (getline(iss, grp, '\t') && getline(iss, clr)) {

                // check if clr is proper rgb value
                istringstream rgb(clr);
                int r, g, b;
                // check if only int given and values are in range 0-255
                if(rgb >> r >> g >> b && r >= 0 && r <= 255 && g >= 0 && g <= 255 && b >= 0 && b <= 255){
                    rgb_color color;
                    color.r = r;
                    color.g = g;
                    color.b = b;

                    grp_clr_map[grp] = color;
                } else {
                    cerr << "Warning: Invalid color (rgb value) given: " << grp << "  " << clr << endl;
                }
            }
        }
        infile.close();

    } else {
        // TODO create colors according to #grpilies
        //  either use color palette (e.g. colorBrewer) or golden ratio
        grp_clr_map.size();
    }

    line = "";
    ifstream plain_nexus(nexus_file);

    // naming output file TODO handle filepath for windeows and unix \\ / and such
    size_t lastSlash = nexus_file.find_last_of("/\\");
    string filename;
    if (lastSlash != std::string::npos) {
        string path = nexus_file.substr(0, lastSlash + 1);
        filename = nexus_file.substr(lastSlash + 1);
        filename = path + "clrd_" + filename;
    } else {
        filename = "clrd_" + nexus_file;
    }
    // colored nexus output file
//    ofstream colored_nexus(filename.c_str());

    if(!plain_nexus.is_open()){
        cerr << "Error opening nexus file to modify with color.\n";
        return;
    }
    /*
    if(!colored_nexus.is_open()){
        cerr << "Error creating colored nexus file " << filename << endl;
        return;
    }

    colored_nexus.close();
    remove(filename.c_str());
    cout << "Created and removed file " << filename << endl;
     */

    // map containing number -> color to color vertices
    unordered_map<int, rgb_color> no_clr_map;

    // read nexus data
    while (getline(plain_nexus, line)) {

        // work with found blocks
        if(translate){

            // String stream to split and read line
            istringstream iss(line);
            string no, taxname; // number taxname of node
            if (getline(iss, no, ' ') && getline(iss, taxname)) {

                // removing quotes and comma
                taxname.erase(remove_if(taxname.begin(), taxname.end(),
                                        [](char c) { return (c == '\'' || c == ','); }), taxname.end());
                // get color via taxname -> group -> color
                rgb_color clr = grp_clr_map[tax_grp_map[taxname]];
                no_clr_map[stoi(no)] = clr;

            } else {
                cerr << "Error reading TRANSLATE block in nexus file. ";
            }

        }
        if(vertices){
            // TODO add color information to line
            // -> VERTICES [number] [vertice] <bg=...> <fg=...>

            // adding color
            //string color = "bg="+bg1+" "+bg2+" "+bg3 +" fg="+fg1+" "+fg2+" "+fg3+",\n";
            //line += color;
        }

        // identify block
        if(line.find("TRANSLATE") != string::npos || line.find("Translate") != string::npos || line.find("translate") != string::npos){// starting point of: [number] [taxname]
            translate = true;

            cout << line << endl;
        }
        if (line.find("VERTICES") != string::npos || line.find("Vertices") != string::npos || line.find("vertices") != string::npos) {// starting point of vertices
            vertices = true;

            cout << line << endl;
        }
        if(line.find(";") != string::npos){// end of any block
            translate = false;
            vertices = false;

            cout << line << endl;
        }

        // colored_nexus << line << endl;
    }

    plain_nexus.close();
    // colored_nexus.close();

    // receive file taxa <-> color
    // read file nexus with already added network via splitstree
    // map color <-> taxa <-> node
    // add color at respective nodes

    // special treatment if nodes would have multiple colors
    //      (draw two dots? other shape? different inner/outer color?)
    // nice-to-have: color for colorblind, or diff shape
}