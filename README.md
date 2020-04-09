# SANS *serif*

### Symmetric Alignment-free phylogeNomic Splits

* Reference-free
* Alignment-free
* Assembled genomes or reads as input
* Phylogenetic splits as output

### Publication

Wittler, R.: [Alignment- and reference-free phylogenomics with colored de Bruijn graphs.](http://drops.dagstuhl.de/opus/volltexte/2019/11032/pdf/LIPIcs-WABI-2019-2.pdf) In: Huber, K. and Gusfield, D. (eds.) Proceedings of WABI 2019. LIPIcs. 143, Schloss Dagstuhl--Leibniz-Zentrum fuer Informatik, Dagstuhl, Germany (2019).



## Table of Contents

* [Requirements](https://gitlab.ub.uni-bielefeld.de/roland.wittler/sans#requirements)
* [Compilation](https://gitlab.ub.uni-bielefeld.de/roland.wittler/sans#compilation)
* [Usage](https://gitlab.ub.uni-bielefeld.de/roland.wittler/sans#usage)
* [FAQ](https://gitlab.ub.uni-bielefeld.de/roland.wittler/sans#faq)
* [Contact](https://gitlab.ub.uni-bielefeld.de/roland.wittler/sans#contact)
* [License](https://gitlab.ub.uni-bielefeld.de/roland.wittler/sans#license)

## Requirements

From the given genomes, a **colored de Bruijn graph** is built to efficiently extract common subsequences. To this end, SANS uses the API of [Bifrost](https://github.com/pmelsted/bifrost). Apart from the requirements of Bifrost (c++ and cmake), there are no further strict dependencies.

* To convert the output into NEXUS format, the provided script requires Phython 3.

* To visualize the splits, we recommend the tool [SplitsTree](http://www.splitstree.org/).


## Compilation

```
cd <sans_directory>
make
```

By default, the installation creates:
* a binary (*SANS*)

You may want to make the binary (*SANS*) accessible via your *PATH* variable.

Please note the installation instructions regsarding the default maximum *k*-mer size of Bifrost in its README.
If your Bifrost libraries have been compiled for 64 bit, change the SANS makefile accordingly (easy to see how).

If during the compilation, the bifrost library files are not found, make sure that the corresponding folder is found as include path by the C++ compiler. You may have to add
`-I/usr/local/include` (with the corresponding folder) to the CFLAGS in the makefile.


## Usage:

```
SANS
```

displays the command line interface:
```

Usage: SANS [PARAMETERS]

  Required arguments:

    -i, --input   	 Input file: list of sequence files, one per line

    -g, --graph   	 Graph file: load a Biforst graph, file name prefix
                  	 (either -i/--input or -g/--graph must be provided)

    -o, --output  	 Output file: list of splits, sorted by weight desc.

  Optional arguments:

    -k, --kmer    	 Length of k-mers (default: 31)

    -t, --top     	 Number of splits (default: all)

    -m, --mean    	 Mean weight function to handle asymmetric splits
                  	 options: arith: arithmetic mean
                  	          geom:  geometric mean (default)
                  	          geom2: geometric mean with pseudo-counts

    -f, --filter  	 Output a greedy maximum weight subset
                  	 options: 1-tree: compatible to a tree
                  	          2-tree: compatible to union of two trees (network)

    -x, --iupac   	 Extended IUPAC alphabet, resolve ambiguous bases
                  	 Specify a number to limit the k-mers per position
                  	 between 1 (no ambiguity) and 4^k (allows NNN...N)

    -v, --verbose 	 Print information messages during execution

    -h, --help    	 Display this help page and quit

```

### Examples

1. **Determine splits from assemblies**
   ```
   SANS -t 4 -k 31 -o sans.splits -r list.txt
   ```
   The colored de Bruijn graph is built with Bifrost using 4 threads (`-t 4`) from the 31-mers (`-k 31`) of those fasta or fastq files listed in *list.txt*  (`-r list.txt`). Splits are determined and written to *sans.splits* (`-o sans.splits`)

   To extract a tree in newick format, use the filter script:
   ```
   scripts/sans2new.py sans.splits > sans_greedytree.new 
   ```

2. **Determine splits from read files**
   ```
   SANS -t 4 -k 31 -o sans.splits -s list.txt
   ```
   The colored de Bruijn graph is built with Bifrost using 4 threads (`-t 4`) from the 31-mers (`-k 31`) of those fasta or fastq files listed in *list.txt*  (`-s list.txt`). By using parameter `-s`, all files are filtered: k-mers occurring exactly once in a file are discarded from the construction.  Splits are determined and written to *sans.splits* (`-o sans.splits`).

3. **Drosophila example data**
   ```
   # go to example directory
   cd <SANS dir>
   cd example_data/drosophila
   
   # download data
   ./download.sh
   
   # run SANS
   cd fa
   SANS -r list.txt -o ../sans.splits -T 130 -t 4 -v
   cd ..
   
   # greedy tree
   ../../scripts/sans2new.py sans.splits -g sans_greedytree.splits > sans_greedytree.new

   # compare to reference
   ../../scripts/newick2sans.py Reference.new > Reference.splits
   ../../scripts/comp.py sans_greedytree.splits Reference.splits fa/list.txt
   ```

4. **Virus example data**
   ```
   # go to example directory
   cd <SANS dir>
   cd example_data/prasinoviruses
      
   # download data
   ./download.sh
   
   # run SANS
   cd fa
   SANS -r list.txt -o ../sans.splits -T 130 -t 4 -v -k11
   cd ..

   # compare to references
   ../../scripts/newick2sans.py Reference_Fig3.new > Reference_Fig3.splits
   ../../scripts/comp.py sans.splits Reference_Fig3.splits fa/list.txt
   ../../scripts/newick2sans.py Reference_Fig4.new > Reference_Fig4.splits
   ../../scripts/comp.py sans.splits Reference_Fig4.splits fa/list.txt
   ```


## FAQ

We recommend to have a look at the [FAQs of Bifrost](https://github.com/pmelsted/bifrost#faq).


## Contact

For any question, feedback or problem, please feel free to file an issue on this Git repository and we will get back to you as soon as possible.

## License

* The hash function library xxHash is BSD licensed (https://github.com/Cyan4973/xxHash)

* The popcount library is BSD licensed (https://github.com/kimwalisch/libpopcnt)

* The libdivide library is zlib licensed (https://github.com/ridiculousfish/libdivide)

* The kseq library is copyrighted by Heng Li and released
  under the MIT license (http://lh3lh3.users.sourceforge.net/kseq.shtml)

* The CRoaring library is Apache 2.0 licensed (https://github.com/RoaringBitmap/CRoaring)

* Bifrost is BSD-2 licensed (https://github.com/pmelsted/bifrost)

* SANS is under [GNU general public license](https://gitlab.ub.uni-bielefeld.de/roland.wittler/sans/blob/master/LICENSE)

