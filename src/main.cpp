#include "main.h"
#include <regex>
// gzstream imports
#include <cstring>
#include "gz/gzstream.h"


/**
 * This is the entry point of the program.
 *
 * @param argc number of cmd args
 * @param argv cmd args
 * @return exit status
 */
int main(int argc, char* argv[]) {
    

    /**
    * [Info]
    * --- Help page ---
    * - Print the help page to console
    * - Describes the program arguments
    */

    if (argc <= 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        cout << endl;
        cout << "SANS serif | version " << SANS_VERSION << endl;
        cout << "Usage: SANS [PARAMETERS]" << endl;
        cout << endl;
        cout << "  Input arguments:" << endl;
        cout << endl;
        cout << "    -i, --input   \t Input file: file of files format" << endl;
        cout << "                  \t Either: one genome per line (space-separated for multifile genomes)" << endl;
        cout << "                  \t Or: kmtricks input format (see https://github.com/tlemane/kmtricks)" << endl;
        cout << endl;
        cout << "    -g, --graph   \t Graph file: load a Bifrost graph, file name prefix" << endl;
#ifndef useBF
        cout << "                  \t (requires compiler flag -DuseBF, please edit makefile)" << endl;
#endif
        cout << endl;
        cout << "    -s, --splits  \t Splits file: load an existing list of splits file" << endl;
        cout << "                  \t (allows to filter -t/-f, other arguments are ignored)" << endl;
        cout << endl;
        cout << "    (either --input and/or --graph, or --splits must be provided)" << endl;
        cout << endl;
        cout << "  Output arguments:" << endl;
        cout << endl;
        cout << "    -o, --output  \t Output TSV file: list of splits, sorted by weight desc." << endl;
        cout << endl;
        cout << "    -N, --newick  \t Output Newick file" << endl;
        cout << "                  \t (only applicable in combination with -f strict or n-tree)" << endl;
        cout << endl;
        cout << "    (at least --output or --newick must be provided, or both)" << endl;
        cout << endl;
        cout << "  Optional arguments:" << endl;
        cout << endl;
        cout << "    -k, --kmer    \t Length of k-mers (default: 31, or 10 for --amino and --code)" << endl;
        cout << endl;
        // cout << "    -w, --window  \t Number of k-mers per minimizer window (default: 1)" << endl;
        // cout << endl;
        cout << "    -t, --top     \t Number of splits in the output list (default: all)." << endl;
        cout << "                  \t Use -t <integer>n to limit relative to number of input files, or" << endl;
        cout << "                  \t use -t <integer> to limit by absolute value." << endl;
        cout << endl;
        cout << "    -m, --mean    \t Mean weight function to handle asymmetric splits" << endl;
        cout << "                  \t options: arith: arithmetic mean" << endl;
        cout << "                  \t          geom:  geometric mean" << endl;
        cout << "                  \t          geom2: geometric mean with pseudo-counts (default)" << endl;
        cout << endl;
        cout << "    -f, --filter  \t Output a greedy maximum weight subset" << endl;
        // cout << "                  \t additional output: (weighted) cleanliness, i.e., ratio of" << endl;
        // cout << "                  \t filtered splits w.r.t. original splits" << endl;
        cout << "                  \t options: strict: compatible to a tree" << endl;
        cout << "                  \t          weakly: weakly compatible network" << endl;
        cout << "                  \t          n-tree: compatible to a union of n trees" << endl;
        cout << "                  \t                  (where n is an arbitrary number, e.g. 2-tree)" << endl;
        cout << endl;
        cout << "    -x, --iupac   \t Extended IUPAC alphabet, resolve ambiguous bases or amino acids" << endl;
        cout << "                  \t Specify a number to limit the k-mers per position between" << endl;
        cout << "                  \t 1 (no ambiguity) and 4^k respectively 22^k (allows NNN...N)" << endl;
        cout << "                  \t Without --iupac respective k-mers are ignored" << endl;
        cout << endl;
        cout << "    -q, --qualify \t Discard k-mers with lower coverage than a threshold" << endl;
        cout << endl;
        cout << "    -n, --norev   \t Do not consider reverse complement k-mers" << endl;
        cout << endl;
        cout << "    -a, --amino   \t Consider amino acids: --input provides amino acid sequences" << endl;
        cout << "                  \t Implies --norev and a default k of 10" << endl;
        cout << endl;
        cout << "    -c, --code   \t Translate DNA: --input provides coding sequences" << endl;
        cout << "                 \t Implies --norev and a default k of 10" << endl;
        cout << "                 \t optional: ID of the genetic code to be used" << endl;
        cout << "                 \t Default: 1" << endl;
        cout << "                 \t Use 11 for Bacterial, Archaeal, and Plant Plastid Code" << endl;
        cout << "                 \t (See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for details.)" << endl;
        cout << endl;
        cout << "    -M, --maxN \t Compare number of input genomes to compile paramter DmaxN" << endl;
        cout << "               \t Add path/to/makefile (default is makefile in current working directory)." << endl;
        cout << endl;
        cout << "    -v, --verbose \t Print information messages during execution" << endl;
        cout << endl;
        cout << "    -h, --help    \t Display this help page and quit" << endl;
        cout << endl;
        cout << "  Contact: sans-service@cebitec.uni-bielefeld.de" << endl;
        cout << "  Evaluation: https://www.surveymonkey.de/r/denbi-service?sc=bigi&tool=sans" << endl;
        cout << endl;
        return 0;
    }


    /**
    * [Meta]
    * --- Defaults ---
    * - Initialise meta variables and set defaults
    */

    string input;    // name of input file
    string graph;    // name of graph file
    string splits;    // name of splits file
    string output;    // name of output file
    string newick;    // name of newick output file // Todo
    string translate; // name of translate file

    uint64_t kmer = 31;    // length of k-mers
    uint64_t window = 1;    // number of k-mers
    uint64_t num = 0;    // number of input files
    uint64_t top = -1;    // number of splits
    bool dyn_top = false; // bind number of splits to num

    auto mean = util::geometric_mean2;    // weight function
    string filter;    // filter function
    uint64_t iupac = 1;    // allow extended iupac characters
    uint64_t quality = 1;    // min. coverage threshold for k-mers
    bool reverse = true;    // consider reverse complement k-mers
    bool verbose = false;    // print messages during execution

    bool amino = false;      // input files are amino acid sequences
    bool shouldTranslate = false;   // translate input files
    bool userKmer = false; // is k-mer default or custom
    bool check_n = false; // compare num (number of input genomes) to maxN (compile parameter DmaxN)
    string path = "./makefile"; // path to makefile
    uint64_t code = 1;


    /**
     * --- Argument parser ---
     * - Parse the command line arguments and update the meta variables accordingly
     */

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) {
            input = argv[++i];    // Input file: list of sequence files, one per line
        }
        else if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--graph") == 0) {
            graph = argv[++i];    // Graph file: load a Bifrost graph, file name prefix
            #ifndef useBF
                cerr << "Error: requires compiler flag -DuseBF" << endl;
                return 1;
            #endif
        }
        else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--splits") == 0) {
            splits = argv[++i];    // Splits file: load an existing list of splits file
        }
        else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            output = argv[++i];    // Output file: list of splits, sorted by weight desc.
        }
        else if (strcmp(argv[i], "-N") == 0 || strcmp(argv[i], "--newick") == 0) {
            newick = argv[++i];    // Output newick file
        }
        else if (strcmp(argv[i], "-k") == 0 || strcmp(argv[i], "--kmer") == 0) {
            kmer = stoi(argv[++i]);    // Length of k-mers (default: 31, 10 for amino acids)
            userKmer = true;
        }
        else if (strcmp(argv[i], "-w") == 0 || strcmp(argv[i], "--window") == 0) {
            window = stoi(argv[++i]);    // Number of k-mers (default: 1)
            if (window > 1) {
                cerr << "Warning: using experimental feature --window" << endl;
            }
        }
        else if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--top") == 0) {
            i++;
            string top_str = argv[i];
            top = stoi(top_str); // Number of splits (default: all)

            if (top_str[top_str.size() - 1] == 'n'){ // Dynamic split num (default: false)
                dyn_top = true;
            }
        }
        else if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--mean") == 0) {
            string arg = argv[++i];    // Mean weight function to handle asymmetric splits
            if (arg == "arith") {
                mean = util::arithmetic_mean;
            }
            else if (arg == "geom") {
                mean = util::geometric_mean;
            }
            else if (arg == "geom2") {
                mean = util::geometric_mean2;
            }
            else {
                cerr << "Error: unknown argument: --mean " << arg << endl;
                return 1;
            }
        }
        else if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--filter") == 0) {
            filter = argv[++i];    // Filter a greedy maximum weight subset
            if (filter == "strict" || filter == "tree") {
                // compatible to a tree
            }
            else if (filter == "weakly") {
                // weakly compatible network
            }
            else if (filter.find("-tree") != -1 && filter.substr(filter.find("-tree")) == "-tree") {
                for (const char &c: filter.substr(0, filter.find("-tree"))){
                    if (!isdigit(c)){
                        cerr << "Error: unexpected argument: --filter " << filter << ". Please specify n (Example usage: --filter 2-tree)" << endl;
                        return 1;
                    }
                }
                stoi(filter.substr(0, filter.find("-tree")));
            }
            else if (filter.find("tree") != -1 && filter.substr(filter.find("tree")) == "tree") {
                for (const char &c: filter.substr(0, filter.find("-tree"))){
                    if (!isdigit(c)){
                        cerr << "Error: unexpected argument: --filter " << filter << ". Please specify n (Example usage: --filter 2-tree)" << endl;
                        return 1;
                    }
                }
                stoi(filter.substr(0, filter.find("tree")));
            }
            else {
                cerr << "Error: unknown argument: --filter " << filter << endl;
                return 1;
            }
        }
        else if (strcmp(argv[i], "-x") == 0 || strcmp(argv[i], "--iupac") == 0) {
            iupac = stoi(argv[++i]);    // Extended IUPAC alphabet, resolve ambiguous bases
        }
        else if (strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "--qualify") == 0) {
            quality = stoi(argv[++i]);    // Discard k-mers below a min. coverage threshold
        }
        else if (strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--norev") == 0) {
            reverse = false;    // Do not consider reverse complement k-mers
        }
        else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;    // Print messages during execution
        }
        else if (strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "--amino") == 0) {
            amino = true;   // Input provides amino acid sequences
        }
        else if (strcmp(argv[i], "-M") == 0 || strcmp(argv[i], "--maxN") == 0) {
            check_n = true; // compare num (number of input genomes) to maxN (compile parameter DmaxN)
            if (i+1 < argc && argv[i+1][0] != '-'){
				path = argv[++i]; // path to makefile
			}
		}
        else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--code") == 0) {
            if (i+1 < argc) {
                string param = argv[++i];
                if (!util::is_number(param)){
                    i--;
                } else {
                    code = stoi(param);     // Number of the genetic code to be used
                }

            }
            shouldTranslate = true;
        }
        else {
            cerr << "Error: unknown argument: " << argv[i] <<  "\t type --help" << endl;
            return 1;
        }
    }

    
    /**
     * --- Version check --- 
     * - Request the current SANS version from gitlab and check if this version is up to date (Requires wget)
     */

    if (verbose){cout << "Checking for updates" << endl;}
    bool version_checked = false;
    if (!system("wget --timeout=1 --tries=1 -qO- https://gitlab.ub.uni-bielefeld.de/gi/sans/raw/master/src/main.h | grep -q SANS_VERSION")){
        version_checked = true;
        if (system("wget --timeout=1 --tries=1 -qO- https://gitlab.ub.uni-bielefeld.de/gi/sans/raw/master/src/main.h | grep -q " SANS_VERSION)) {
        cout << "NEW VERSION AVAILABLE: https://gitlab.ub.uni-bielefeld.de/gi/sans" << endl;
        }
        else if(verbose){cout << "Version up to date" << endl;}
    }
    if (!version_checked && verbose) {cout << "Could not fetch version information" << endl;}


    /**
     * --- Restriction check ---
     * - Check if the given argument configuration does violate any run restrictions
     */ 
    if (!userKmer) {
        if (amino || shouldTranslate) {
            kmer = 10;
        } else {
            kmer = 31;
        }
    }

    if (input.empty() && graph.empty() && splits.empty()) {
        cerr << "Error: missing argument: --input <file_name> or --graph <file_name>" << endl;
        return 1;
    }
    if (!input.empty() && !graph.empty() && !splits.empty()) {
        cerr << "Error: too many input arguments: --input, --graph, and --splits" << endl;
        return 1;
    }
    if (!graph.empty() && !splits.empty()) { // ---- Why not?
        cerr << "Error: too many input arguments: --graph and --splits" << endl;
        return 1;
    }

    if (input.empty() && amino) {
        cerr << "Error: missing argument: --input <file_name> for option --amino" << endl;
        return 1;
    }
    if (!splits.empty() && amino) {
        cerr << "Error: too many input arguments: --splits and --amino" << endl;
        return 1;
    }
    if (!graph.empty() && amino) {
        cerr << "Error: too many input arguments: --graph and --amino" << endl;
        return 1;
    }

    if (output.empty() && newick.empty()) {
        cerr << "Error: missing argument: --output <file_name> or --newick <file_name>" << endl;
        return 1;
    }
    if (kmer > maxK && splits.empty()) {
        cerr << "Error: k-mer length exceeds -DmaxK=" << maxK << endl;
        return 1;
    }
    if (!newick.empty() && filter != "strict" && filter.find("tree") == -1) {
        cerr << "Error: Newick output only applicable in combination with -f strict or n-tree" << endl;
        return 1;
    }

    if (amino && quality > 1) {
        cerr << "Error: using --qualify with --amino is (currently) not supported" << endl;
        cerr << "       Please send us a message if you have this special use case" << endl;
        return 1;
    }
    if (!graph.empty()) {
        if (quality > 1) {
            cerr << "Warning: input from graph with --qualify can produce unexpected results" << endl;
        }
        else if (window > 1) {
            cerr << "Warning: input from graph with --window can produce unexpected results" << endl;
        }
    } else if (quality > 1 && window > 1) {
        cerr << "Warning: using --qualify with --window could produce unexpected results" << endl;
    }

    if (!input.empty() && !splits.empty()) {
        cerr << "Note: two input arguments --input and --splits were provided" << endl;
        cerr << "      --input is used for lookup only, no additional splits are inferred" << endl;
    }
    if (input.empty() && !splits.empty() && !newick.empty()) {
        cerr << "Note: Newick output from a list of splits, some taxa could be missing" << endl;
        cerr << "      --input can be used to provide the original list of taxa" << endl;
    }

    // check if we need to init translation
    if (shouldTranslate) {
        amino = true;
        if (!translator::init(code)) {
            cerr << "Error: No translation data found" << translate << endl;
        }
    }

    /**
     * [Indexing]
     * --- Indexing inputs ---
     * - Collect all target file names and locations
     * - Check if the given files exist
     */ 

    // determine the folder the list is contained in
    string folder="";
	uint64_t found=input.find_last_of("/\\");
	if (found!=string::npos){
	  folder=input.substr(0,found+1);
	}

    // parse the list of input sequence files
    hash_map<string, uint64_t> name_table; // the name to color map
    vector<string> denom_names; // storing the representative name per color
    vector<vector<string>> gen_files; // genome file collection

    if (!input.empty()) {
        // check the input file 
        ifstream file(input);
        if (!file.good()) {
            cerr << "Error: could not read input file: " << input << endl;
            return 1;
        }


        // parse the list of input sequence files
        string line; // the iterated input line
        string file_name; // the current file name
        bool is_first; // indicating the first filename of a line (For file list)
        bool has_files; // indicating if a line contains filenames

        getline(file, line);
        // check the file format
        std::smatch matches;
        std::regex_search(line, matches, std::regex("(:)"));
        bool is_kmt = !matches.empty();

        // parse file of files
        while(true){
            vector<string> target_files; // container of the current target files
            if (is_kmt){ // parse kmt format
                // ensure the terminal signs " !" exists.
                if (line.find_first_of('!') != line.npos){line = line.substr(0, line.find_first_of('!') + 1);} // cut off unused tail
                else if (line.back() == ' '){line += '!';} // append terminal sign if missing
                else {line += " !";} // append both terminal signs if missing

                string denom = line.substr(0, line.find_first_of(" ")); // get the dataset-id
                denom_names.push_back(denom); // add id to denominators
                num ++;
                line = line.substr(line.find_first_of(":") + 2, line.npos); // cut off the dataset-id

                std::smatch matches; // Match files
                while (std::regex_search(line, matches, std::regex("[ ; ]|[ !]"))){
                    file_name = matches.prefix().str(); // get filename from match
                    line=matches.suffix().str(); // update the line
                    if (file_name.length() == 0){continue;} // skip empty file name
                    else {target_files.push_back(file_name); name_table[file_name] = num;} // add the file name to target files and name table
                    }
            }

            else{ // parse file list format
                is_first = true;
                has_files = false;
                string file_name = "";
                size_t it = 0;
                size_t line_length = line.length();
                for (auto x: line){ // iterate the line
                    it ++;
                    if (x == ' ' | it == line_length){ // checkout the file name if a space occurs or the line ends
                        if (it == line_length){file_name += x;} // add the last character to the last file name
                        if (file_name.length() == 0){file_name = ""; continue;} // skip continuous spaces
                        if (is_first){ // use first file name as denom name
                            has_files = true;
                            denom_names.push_back(file_name); // set denom name
                            is_first = false;
                        }
                        target_files.push_back(file_name); // add the file_name to the genome file vector
                        name_table[file_name] = num; // add the file tp the name_table
                        file_name = "";
                    }
                    else{file_name += x;}
                }
                if (has_files) {num++;}
            }


            // check files
	    if (splits.empty()){
            	for(string file_name: target_files){
                    ifstream file_stream = ifstream(folder+file_name);
                    if (!file_stream.good()) { // catch unreadable file
                        cout << "\33[2K\r" << "\u001b[31m" << "(ERR)" << " Could not read file " <<  "<" << folder+file_name << ">" << "\u001b[0m" << endl;
                    	file_stream.close();
                    	return 1;
                    }
                    else{ file_stream.close();}	
            	}
	    }
            gen_files.push_back(target_files); // add the files of the current genome to the genome collection
            if (!getline(file, line)) {break;}
        }
    }
    int denom_file_count = denom_names.size(); 


    /**
     * --- Indexing CDBG input --- 
     * - Collect all target sequence names from the Bifrost CDBG
     */ 

#ifdef useBF
    // load an existing Bifrost graph
    ColoredCDBG<> cdbg(kmer);
    if (!graph.empty()) {
        if (cdbg.read(graph + ".gfa", graph + ".bfg_colors", 1, verbose)) { // Allow parallel reading with new t parameter.

        } else {
            cerr << "Error: could not load Bifrost graph" << endl;
            return 1;
        }

        if (kmer != cdbg.getK()) {
            kmer = cdbg.getK();
            if (kmer > maxK) {
                cerr << "Error: k-mer length exceeds -DmaxK=" << maxK << endl;
                return 1;
            }
            cerr << "Warning: setting k-mer length to " << kmer << endl;
        }

        vector<string> cdbg_names = cdbg.getColorNames(); // color names of the cdbg compacted genomes.
        for (auto it=0; it != cdbg_names.size(); ++it){ // iterate the cdbg names and transcribe them to the name table
            if (name_table.find(cdbg_names[it]) == name_table.end()){
                            name_table[cdbg_names[it]] = num++;
                            denom_names.push_back(cdbg_names[it]);
			    vector<string> dummy;	
			    gen_files.push_back(dummy);
	        }
            else
            {
                cout << "Warning: " << cdbg_names[it] << " exists in input and graph. It is treated as one sequence" << endl;
            }
        }

        if (verbose) {
            cout << endl;
	}
    }

#endif


    /**
     * --- Post indexing check ---
     * - Update and check validity of input dependent meta variables
     */ 

    // check if the number of genomes is reasonably close the maximal storable color set
    if (check_n) {
       util::check_n(num,path);
	}

    // check if the number of gernomes exceeds the maximal storable color set
    if (num > maxN) {
        cerr << "Error: number of input genomes ("<<num<<") exceeds -DmaxN=" << maxN << endl;
        return 1;
    }

    // Set dynamic top by filenum
    if (dyn_top){
        top = top * num;
    }


    /**
     * [Input processing]
     * - Transcribe each given input to the graph
     */ 

    chrono::high_resolution_clock::time_point begin = chrono::high_resolution_clock::now();    // time measurement
    graph::init(top, quality, amino); // initialize the top list size, coverage threshold, and allowed characters


    /**
     *  --- Split processing ---
     *  - Transcibe all splits from the input split file
     */

    // iterate splits and add them to the toplist
    if (!splits.empty()) {
    ifstream file(splits);
    if (!file.good()) { // check if the target splits file exists
        cerr << "Error: could not read splits file: " << splits << endl;
        return 1;
    }
    string line;
    while (getline(file, line)) { // Iterate each split
        uint64_t curr = line.find('\t');
        double weight = stod(line.substr(0, curr));
        uint64_t next = curr + 1;

        color_t color = 0;
        do {
            curr = line.find('\t', next);
            string name = line.substr(next, curr-next);
            if (name_table.find(name) == name_table.end()) { // check if the splits genome names are already indexed
                vector<string> file_vec;
                name_table[name] = num++;
                denom_names.push_back(name);
            }
            color::set(color, name_table[name]);
            next = curr + 1;
        } while (curr != string::npos);

        graph::add_split(weight, color); // add the split to the graph
    }
    file.close();
    }

    /**
     * --- Sequence processing ---
     * - Translate all given sequence k-mers
     * - Transcribe all given sequence k-mers to the graph
     */ 
    
    kmer::init(kmer);      // initialize the k-mer length
    kmerAmino::init(kmer); // initialize the k-mer length
    color::init(num);    // initialize the color number

    if (!input.empty() && splits.empty()) {
        if (verbose) {
            cout << "SANS::main(): Reading input files..." << flush;
        }
        string sequence;    // read in the sequence files and extract the k-mers


        for (uint64_t i = 0; i < gen_files.size(); ++i) {
            vector<string> target_files = gen_files[i]; // the filenames corresponding to the target	    
            for (string file_name: target_files){
		
	 	        char c_name[(folder + file_name).length()]; // Create char array for c compatibilty
		        strcpy(c_name, (folder + file_name).c_str()); // Transcire to char array

                igzstream file(c_name, ios::in);    // input file stream
                if (verbose) {
                    cout << "\33[2K\r" << folder+file_name<< " (" << i+1 << "/" << denom_file_count << ")" << endl;    // print progress
                }
                count::deleteCount();

                string appendixChars; 
                string line;    // read the file line by line
                while (getline(file, line)) {
                    if (line.length() > 0) {
                        if (line[0] == '>' || line[0] == '@') {    // FASTA & FASTQ header -> process
                            if (window > 1) {
                                iupac > 1 ? graph::add_minimizers(sequence, i, reverse, window, iupac)
                                        : graph::add_minimizers(sequence, i, reverse, window);
                            } else {
                                iupac > 1 ? graph::add_kmers(sequence, i, reverse, iupac)
                                        : graph::add_kmers(sequence, i, reverse);
                            }

                            sequence.clear();

                            if (verbose) {
                                cout << "\33[2K\r" << line << flush << endl;    // print progress
                            }
                        }
                        else if (line[0] == '+') {    // FASTQ quality values -> ignore
                            getline(file, line);
                        }
                        else {
                            transform(line.begin(), line.end(), line.begin(), ::toupper);
                            string newLine = line;
                            if (shouldTranslate) {
                                if (appendixChars.length() >0 ) {
                                    newLine= appendixChars + newLine;
                                    appendixChars = "";
                                }
                                auto toManyChars = line.length() % 3;
                                if (toManyChars > 0) {
                                    appendixChars = newLine.substr(line.length() - toManyChars, toManyChars);
                                    newLine = newLine.substr(0, line.length() - toManyChars);
                                }

                                newLine = translator::translate(newLine);
                            }
                            sequence += newLine;    // FASTA & FASTQ sequence -> read
                        }
                    }
                }
                if (verbose && count::getCount() > 0) {
                    cerr << count::getCount()<< " triplets could not be translated."<< endl;
                }
                if (window > 1) {
                    iupac > 1 ? graph::add_minimizers(sequence, i, reverse, window, iupac)
                            : graph::add_minimizers(sequence, i, reverse, window);
                } else {
                    iupac > 1 ? graph::add_kmers(sequence, i, reverse, iupac)
                            : graph::add_kmers(sequence, i, reverse);
                }
                sequence.clear();

                
                if (verbose) {
                    cout << "\33[2K\r" << flush;
                }
                file.close();
            }
            graph::clear_thread();
        }
    }


    /**
     * --- Bifrost CDBG processing ---
     * - Iterate all colored k-mers from a CDBG
     * - Compute the splits created by the CDBG k-mers given the graphs colore k-mer collection
     * (Has to be executed after sequence processing)
     */ 

double min_value = numeric_limits<double>::min(); // Current minimal weight represented in the top list
#ifdef useBF
    if (!graph.empty()) {
        if (verbose) {
            cout << "SANS::main(): Processing unitigs..." << flush;
        }
        uint64_t cur = 0, progress;
        uint64_t max = cdbg.size();

        for (auto& unitig : cdbg) {
            if (verbose) {
                cout << "\33[2K\r" << "Processed " << cur << " unitigs (" << 100*cur/max << "%) " << flush;
            }   cur++;

            auto num_kmers = unitig.size - kmer::k + 1; // the number of kmers in this unitig
            auto uc_kmers = new UnitigColors[num_kmers]; // storage for unitig colored kmers
            auto unitig_map = UnitigMapBase(0, 1, kmer::k, true);

            auto sequence = unitig.mappedSequenceToString(); // the mapped unitig sequence
            auto *colors = unitig.getData()->getUnitigColors(unitig); // the k-mer-position-per-color iterator of this unitig
            auto it = colors->begin(unitig);
            auto end = colors->end();

            for (; it != end; ++it) { // iterate the unitig and collect the colors and corresponding k-mer starts
                uc_kmers[it.getKmerPosition()].add(unitig_map, it.getColorID());
            }
            
            for (unsigned int i = 0; i != num_kmers; ++i){ // iterate the k-mers
                string kmer_sequence = sequence.substr(i, kmer::k); // the k-mer sequence
                color_t color = 0;
                for (auto uc_it=uc_kmers[i].begin(unitig_map); uc_it != uc_kmers[i].end(); ++uc_it){
                    color::set(color, name_table[cdbg.getColorName(uc_it.getColorID())]); // set the k-mer color
		}
                min_value = graph::add_cdbg_colored_kmer(mean, kmer_sequence, color, min_value);
	   }
        }
    }
#endif


    /**
     * [Output]
     * - Compute the weighted splits of the k-mers that are left in the graph
     * - Apply the target filter method
     * - Write to output
     */ 

    // function to map color position to file name
    std::function<string(const uint64_t&)> map=[=](uint64_t i) {
        if (i < denom_names.size()) return denom_names[i];
        cerr << "Error: color bit does not correspond to color name" << endl;
        exit(EXIT_FAILURE);
    };

    if (verbose) {
        cout << "Processing splits..." << flush;
    }
    graph::add_weights(mean, min_value, verbose);  // accumulate split weights

    if (verbose) {
        cout << "\33[2K\r" << "Filtering splits..." << flush;
    }

    //cleanliness cleanliness;
    if (!filter.empty()) {    // apply filter
        for (auto& split : graph::split_list) {
            //cleanliness.addWeightStateBefore(split.first, split.second);
        }

        if (filter == "strict" || filter == "tree") {
            if (!newick.empty()) {
                ofstream file(newick);    // output file stream
                ostream stream(file.rdbuf());
                stream << graph::filter_strict(map, verbose);    // filter and output
                file.close();
            } else {
               graph::filter_strict(verbose);
            }
        }
        else if (filter == "weakly") {
           graph::filter_weakly(verbose);
        }
        else if (filter.find("tree") != -1 && filter.substr(filter.find("tree")) == "tree") {
            if (!newick.empty()) {
                ofstream file(newick);    // output file stream
                ostream stream(file.rdbuf());
                auto ot = graph::filter_n_tree(stoi(filter.substr(0, filter.find("tree"))), map, verbose);
                stream <<  ot;
                file.close();
            } else {
               graph::filter_n_tree(stoi(filter.substr(0, filter.find("tree"))), verbose);
            }
        }
    }


    /**
     * --- Write to output ---
     */ 

    if (verbose) {
        cout << "\33[2K\r" << "Please wait..." << flush << endl;
    }
    ofstream file(output);    // output file stream
    ostream stream(file.rdbuf());

    uint64_t pos = 0;
    //cleanliness.setFilteredCount(graph::split_list.size());
    for (auto& split : graph::split_list) {
        double weight = split.first;
       // cleanliness.setSmallestWeight(weight, split.second);
        stream << weight;    // weight of the split
        for (uint64_t i = 0; i < num; ++i) {
            if (color::test(split.second, pos)) {
                if (i < denom_names.size())
                    stream << '\t' << denom_names[i];    // name of the file
            }
            split.second >>= 01u;
        }
        stream << endl;
    }

    //cleanliness.calculateWeightBeforeCounter();

    file.close();

    chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();    // time measurement

    if (verbose) {
        if (!filter.empty()) {
          // cleanliness.reportCleanliness();
        }
        cout << " Done!" << flush << endl;    // print progress and time
        cout << " (" << util::format_time(end - begin) << ")" << endl;
    }
    return 0;
}
