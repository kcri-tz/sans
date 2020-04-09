#include "main.h"

/**
 * This is the entry point of the program.
 *
 * @param argc number of cmd args
 * @param argv cmd args
 * @return exit status
 */
int main(int argc, char* argv[]) {

    // print a help message describing the program arguments
    if (argc <= 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        cout << endl;
        cout << "Usage: SANS [PARAMETERS]" << endl;
        cout << endl;
        cout << "  Required arguments:" << endl;
        cout << endl;
        cout << "    -i, --input   \t Input file: list of sequence files, one per line" << endl;
        cout << endl;
        cout << "    -g, --graph   \t Graph file: load a Biforst graph, file name prefix" << endl;
        cout << "                  \t (either -i/--input or -g/--graph must be provided)" << endl;
        cout << endl;
        cout << "    -o, --output  \t Output file: list of splits, sorted by weight desc." << endl;
        cout << endl;
        cout << "  Optional arguments:" << endl;
        cout << endl;
        cout << "    -k, --kmer    \t Length of k-mers (default: 31)" << endl;
        cout << endl;
        cout << "    -t, --top     \t Number of splits (default: all)" << endl;
        cout << endl;
        cout << "    -m, --mean    \t Mean weight function to handle asymmetric splits" << endl;
        cout << "                  \t options: arith: arithmetic mean" << endl;
        cout << "                  \t          geom:  geometric mean (default)" << endl;
        cout << "                  \t          geom2: geometric mean with pseudo-counts" << endl;
        cout << endl;
        cout << "    -f, --filter  \t Output a greedy maximum weight subset" << endl;
        cout << "                  \t options: 1-tree: compatible to a tree" << endl;
        cout << "                  \t          2-tree: compatible to union of two trees (network)" << endl;
        cout << endl;
        cout << "    -x, --iupac   \t Extended IUPAC alphabet, resolve ambiguous bases" << endl;
        cout << "                  \t Specify a number to limit the k-mers per position" << endl;
        cout << "                  \t between 1 (no ambiguity) and 4^k (allows NNN...N)" << endl;
        cout << endl;
        cout << "    -v, --verbose \t Print information messages during execution" << endl;
        cout << endl;
        cout << "    -h, --help    \t Display this help page and quit" << endl;
        cout << endl;
        return 0;
    }

    string input;    // name of input file
    string graph;    // name of graph file
    string output;    // name of output file
    uint64_t num = 0;    // number of input files

    uint64_t kmer = 31;    // length of k-mers
    uint64_t top = 0;    // number of splits
    auto mean = util::geometric_mean;    // weight function
    auto filter = graph::filter_none;    // filter function
    uint64_t iupac = 0;    // allow extended iupac characters
    bool verbose = false;    // print messages during execution

    // parse the command line arguments and update the variables above
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) {
            input = argv[++i];    // Input file: list of sequence files, one per line
        }
        else if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--graph") == 0) {
            graph = argv[++i];    // Graph file: load a Biforst graph, file name prefix
            #ifndef useBF
                cerr << "Error: requires compiler flag -DuseBF" << endl;
                return 1;
            #endif
        }
        else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            output = argv[++i];    // Output file: list of splits, sorted by weight desc.
        }
        else if (strcmp(argv[i], "-k") == 0 || strcmp(argv[i], "--kmer") == 0) {
            kmer = stoi(argv[++i]);    // Length of k-mers (default: 31)
        }
        else if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--top") == 0) {
            top = stoi(argv[++i]);    // Number of splits (default: all)
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
            string arg = argv[++i];    // Filter a maximum weight n-tree compatible subset
            if (arg == "none") {
                filter = graph::filter_none;
            }
            else if (arg == "1-tree") {
                filter = graph::filter_tree1;
            }
            else if (arg == "2-tree") {
                filter = graph::filter_tree2;
            }
            else {
                cerr << "Error: unknown argument: --filter " << arg << endl;
                return 1;
            }
        }
        else if (strcmp(argv[i], "-x") == 0 || strcmp(argv[i], "--iupac") == 0) {
            iupac = stoi(argv[++i]);    // Extended IUPAC alphabet, resolve ambiguous bases
        }
        else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;    // Print messages during execution
        }
        else {
            cerr << "Error: unknown argument: type --help" << endl;
            return 1;
        }
    }

    if (input.empty() && graph.empty()) {
        cerr << "Error: missing argument: --input <file_name>" << endl;
        return 1;
    }
    if (output.empty()) {
        cerr << "Error: missing argument: --output <file_name>" << endl;
        return 1;
    }
    if (kmer > maxK) {
        cerr << "Error: k-mer length exceeds -DmaxK=" << maxK << endl;
        return 1;
    }

    // parse the list of input sequence files
    vector<string> files;
    if (!input.empty()) {
        ifstream file(input);
        string line;
        while (getline(file, line)) {
            files.emplace_back(line);
            num++;
        }
        file.close();

        if (num > maxN) {
            cerr << "Error: number of files exceeds -DmaxN=" << maxN << endl;
            return 1;
        }
    }

#ifdef useBF
    // load an existing Bifrost graph
    ColoredCDBG<> cdbg(kmer);
    if (!graph.empty()) {
        if (cdbg.read(graph + ".gfa", graph + ".bfg_colors", 1, verbose)) {
            num += cdbg.getNbColors();
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
        if (num > maxN) {
            cerr << "Error: number of colors exceeds -DmaxN=" << maxN << endl;
            return 1;
        }
        if (verbose) {
            cout << endl;
        }
    }
#endif

    kmer::init(kmer);    // initialize the k-mer length
    color::init(num);    // initialize the color number
    graph::init(top);    // initialize the toplist size

    chrono::high_resolution_clock::time_point begin = chrono::high_resolution_clock::now();    // time measurement

    if (!input.empty()) {
        if (verbose) {
            cout << "SANS::main(): Reading input files..." << flush;
        }
        string sequence;    // read in the sequence files and extract the k-mers

        for (uint64_t i = 0; i < files.size(); ++i) {
            if (verbose) {
                cout << "\33[2K\r" << files[i] << " (" << i+1 << "/" << files.size() << ")" << endl;    // print progress
            }
            ifstream file(files[i]);    // input file stream

            string line;    // read the file line by line
            while (getline(file, line)) {
                if (line.length() > 0) {
                    if (line[0] == '>' || line[0] == '@') {    // FASTA & FASTQ header -> process
                        transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
                        iupac > 0 ? graph::add_kmers_iupac(sequence, i, iupac)
                                  : graph::add_kmers(sequence, i);
                        sequence.clear();

                        if (verbose) {
                            cout << "\33[2K\r" << line << flush;    // print progress
                        }
                    }
                    else if (line[0] == '+') {    // FASTQ quality values -> ignore
                        getline(file, line);
                    }
                    else {
                        sequence += line;    // FASTA & FASTQ sequence -> read
                    }
                }
            }
            transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
            iupac > 0 ? graph::add_kmers_iupac(sequence, i, iupac)
                      : graph::add_kmers(sequence, i);
            sequence.clear();

            if (verbose) {
                cout << "\33[2K\r" << flush;
            }
            file.close();
        }
    }

#ifdef useBF
    if (!graph.empty()) {
        if (verbose) {
            cout << "SANS::main(): Processing unitigs..." << flush;
        }
        uint64_t cur = 0;
        uint64_t max = cdbg.size();

        for (auto& unitig : cdbg) {
            if (verbose) {
                cout << "\33[2K\r" << "Processed " << cur << " unitigs (" << 100*cur/max << "%) " << flush;
            }
            cur++;

            auto sequence = unitig.mappedSequenceToString();
            auto *colors = unitig.getData()->getUnitigColors(unitig);
            auto it = colors->begin(unitig);
            auto end = colors->end();

            for (; it != end; ++it) {
                auto seq = sequence.substr(it.getKmerPosition(), kmer);
                auto col = files.size() + it.getColorID();
                graph::add_kmers(seq, col);
            }
        }
        if (verbose) {
            cout << "\33[2K\r" << "Processed " << max << " unitigs (100%)" << endl;
        }
    }
#endif

    if (verbose) {
        cout << "Please wait." << flush;
    }
    graph::add_weights(mean);    // accumulate split weights

    if (verbose) {
        cout << "." << flush;
    }
    filter();    // apply n-tree filter

    if (verbose) {
        cout << "." << flush;
    }
    ofstream file(output);    // output file stream
    ostream stream(file.rdbuf());

    uint64_t pos = 0;
    for (auto& split : graph::split_list) {
        stream << split.first;    // weight of the split
        for (uint64_t i = 0; i < num; ++i) {
            if (color::test(split.second, pos)) {
                if (i < files.size())
                    stream << "\t" << files[i];    // name of the file
                #ifdef useBF
                else
                    stream << "\t" << cdbg.getColorName(i-files.size());
                #endif
            }
            split.second >>= 01u;
        }
        stream << endl;
    }
    file.close();

    chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();    // time mesaurement

    if (verbose) {
        cout << " Done!" << flush;    // print progress and time
        cout << " (" << util::format_time(end - begin) << ")" << endl;
    }

    return 0;
}
