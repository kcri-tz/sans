#include "nexus_color.h"

struct rgb_color {
    int r, g, b;

    bool is_white(){
        return (r == 0) && (g == 0) && (b == 0);
    }
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
 * Also creates another map with the given values as keys for further usage.
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

// Using the golden ratio to create distinct hsv colors and then convert them to rgb
vector<rgb_color> create_colors(int n) { // !not my function!
    vector<rgb_color> colors;

    // Golden ratio to ensure visual distinctiveness
    const double goldenRatioConjugate = 0.618033988749895;
    double hue = 0.0;

    for (int i = 0; i < n; ++i) {
        // Convert HSV to RGB
        double saturation = 0.7;
        double value = 0.9;

        int hi = static_cast<int>(std::floor(hue * 6.0));
        double f = hue * 6.0 - hi;
        double p = value * (1.0 - saturation);
        double q = value * (1.0 - f * saturation);
        double t = value * (1.0 - (1.0 - f) * saturation);

        int r, g, b;
        switch (hi % 6) {
            case 0: r = value * 255; g = t * 255; b = p * 255; break;
            case 1: r = q * 255; g = value * 255; b = p * 255; break;
            case 2: r = p * 255; g = value * 255; b = t * 255; break;
            case 3: r = p * 255; g = q * 255; b = value * 255; break;
            case 4: r = t * 255; g = p * 255; b = value * 255; break;
            case 5: r = value * 255; g = p * 255; b = q * 255; break;
        }

        rgb_color color;
        color.r = r;
        color.g = g;
        color.b = b;
        colors.push_back(color);

        // Increment hue using golden ratio
        hue += goldenRatioConjugate;
        hue = fmod(hue, 1.0);
    }

    return colors;
}


void nexus_color::open_in_splitstree(const string& nexus_file, const string& pdf, bool verbose, bool update, const string splitstree_path){
    // = "../splitstree4/SplitsTree"
    /*
     * begin SplitsTree;
EXECUTE FILE=example_data/pra_test.nex
EXPORTGRAPHICS format=PDF file=test.pdf REPLACE=yes
QUIT
end;

     */
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
        if(update) temp_file << "UPDATE\nSAVE FILE=" << nexus_file << " REPLACE=YES\n";
        if(!pdf.empty()) temp_file << "EXPORTGRAPHICS format=PDF file=" << pdf << " REPLACE=yes\n";
        temp_file << "QUIT\nend;";
        temp_file.close();

        // Running splitstree
        if(verbose) cout  << "Running SplitsTree\n";
        /// .../SplitsTree -g -S -c run_splitstree
        string command = splitstree_path + " -g -S -c " + temp;
        const char* splitstree_command = command.c_str();
        int result = system(splitstree_command); // executing command
        if(result != 0){
            cerr << "Error while running Splitstree " << splitstree_path << endl;
        }
        remove(temp);
    } else {
        cerr << "Unable to create temp command file for splitstree: " << temp << "\nTherefor unable to add network to nexus file" << endl;
    }
}


void nexus_color::color_nexus(const string& nexus_file, const string& tax_grp_file, const string& grp_clr_file){

    bool translate = false;
    bool vertices = false;
    //vector<string> lines;
    string line;

    unordered_map<string, string> tax_grp_map; // map containing taxname -> group
    unordered_map<string, rgb_color> grp_clr_map; // map conatining group -> color
    unordered_map<int, rgb_color> no_clr_map; // map containing number -> color to color vertices

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

                istringstream rgb(clr);
                int r, g, b;
                // check if clr is proper rgb value: only ints and in range 0-255
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

    } else { // create colors if they are not given

        int no_grps = grp_clr_map.size(); // number of groups
        vector<rgb_color> colors = create_colors(no_grps);

        int i = 0;
        for(auto& vertex : grp_clr_map){
            vertex.second = colors[i];
            ++i;
        }
    }

    line = "";
    ifstream plain_nexus(nexus_file);

    // naming output file TODO handle filepath for windeows and unix \\ / and such
    // TODO necessary??
    size_t lastSlash = nexus_file.find_last_of("/\\");
    string filename;
    if (lastSlash != string::npos) {
        string path = nexus_file.substr(0, lastSlash + 1);
        filename = nexus_file.substr(lastSlash + 1);
        filename = path + "clrd_" + filename;
    } else {
        filename = "clrd_" + nexus_file;
    }

    // colored nexus output file
    ofstream colored_nexus(filename.c_str());

    if(!plain_nexus.is_open()){
        cerr << "Error opening nexus file to modify with color.\n";
        return;
    }
    if(!colored_nexus.is_open()){
        cerr << "Error creating colored nexus file " << filename << endl;
        return;
    }


    // read nexus data
    while (getline(plain_nexus, line)) {

        // end of any block
        if(line.find(";") != string::npos){
            translate = false;
            vertices = false;
        }

        // work with found blocks
        if(translate){
            istringstream iss(line);
            string no, taxname; // number, taxname of node
            if (getline(iss, no, ' ') && getline(iss, taxname)) {

                // removing quotes and comma
                taxname.erase(remove_if(taxname.begin(), taxname.end(),
                                        [](char c) { return (c == '\'' || c == ','); }), taxname.end());
                // get color via taxname -> group -> color
                rgb_color clr = grp_clr_map[tax_grp_map[taxname]];
                no_clr_map[stoi(no)] = clr;

            } else {
                cerr << "Error reading TRANSLATE block in nexus file.\n";
                cerr << "Problematic line: " << line << endl;
            }

        }
        if(vertices) {
            // variables to save info contained for the vertex
            istringstream iss(line);
            int vertex_no;
            double x, y;
            char c_w, c_h, c_s, s, eq;
            int w, h;
            string rest;
            int size = 8; // TODO anpassen ?
            // splitting line into its components: [no.] [x] [y] w=[..] h=[..] s=n,
            if (iss >> vertex_no >> x >> y >> c_w >> eq >> w >> c_h >> eq >> h >> c_s >> eq >> s >> rest) {

                if (no_clr_map.find(vertex_no) != no_clr_map.end()){ // only color vertex if leaf

                    if(rest == ","){ // if last char is ',' add color
                        rgb_color clr = no_clr_map[vertex_no]; // get defined color for vertex
                        int bg_r, bg_g, bg_b;
                        // add color information to line
                        ostringstream oss;
                        if(clr.is_white()){ // if color is white, turn outline black
                            bg_r = bg_g = bg_b = 255;
                        } else {
                            bg_b = clr.b;
                            bg_g = clr.g;
                            bg_r = clr.r;
                        }
                        oss <<vertex_no<<" "<<x<<" "<<y<<" w="<<size<<" h="<<size<<" fg="<<clr.r<<" "<<clr.g<<" "<< clr.b<<" bg="<<bg_r<<" "<<bg_g<<" "<<bg_b<<",";
                        line = oss.str();

                    } else { /* if already colored do not change */ }

                } else {
                    // TODO check what it does for lines that should not be colored
                }
            } else {
                // TODO
                cout << "Problems reading vertex: " << line << endl;
            }
        }

        // identify block
        if(line.find("TRANSLATE") != string::npos || line.find("Translate") != string::npos || line.find("translate") != string::npos){// starting point of: [number] [taxname]
            translate = true;
        }
        if (line.find("VERTICES") != string::npos || line.find("Vertices") != string::npos || line.find("vertices") != string::npos) {// starting point of vertices
            vertices = true;
        }

        colored_nexus << line << endl;
    }

    plain_nexus.close();
    colored_nexus.close();

    std::rename(filename.c_str(), nexus_file.c_str());
    // special treatment if nodes would have multiple colors
    //      (draw two dots? other shape? different inner/outer color?)
    // nice-to-have: color for colorblind, or diff shape
}