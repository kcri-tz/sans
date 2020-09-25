#include "graph.h"

/**
 * This is the size of the top list.
 */
uint64_t graph::t;

/**
 * This is a hash table mapping k-mers to colors [O(1)].
 */
unordered_map<kmer_t, color_t> graph::kmer_table;

/**
 * This is a hash table mapping colors to weights [O(1)].
 */
unordered_map<color_t, array<uint32_t,2>> graph::color_table;

/**
 * This is an ordered tree collecting the splits [O(log n)].
 */
multimap<double, color_t, greater<>> graph::split_list;



/**
 * Initializes a new set struct.
 * @param taxa color_t coding all taxa beneath this node
 * @param subsets list of subsets
 */
struct set* newSet(color_t taxa, double weight, vector<set *> subsets) { 
	// declare and allocate new node  
	struct set* set = new struct set(); 
	set->taxa=taxa;
	set->weight=weight;
	set->subsets = subsets;
	return(set); 
} 


/**
 * This function initializes the top list size.
 *
 * @param t top list size
 */
void graph::init(uint64_t& top_size) {
    t = top_size;
}

/**
 * This function extracts k-mers from a sequence and adds them to the hash table.
 *
 * @param str dna sequence
 * @param color color flag
 * @param reverse merge complements
 */
void graph::add_kmers(string& str, uint64_t& color, bool& reverse) {

    uint64_t pos;    // current position in the string, from 0 to length
    kmer_t kmer;    // create a new empty bit sequence for the k-mer
    kmer_t rcmer;    // create a bit sequence for the reverse complement

next_kmer:
    if (str.length() < kmer::k) {
        return;    // not enough characters for a k-mer, abort
    }
    pos = 0;

    for (; pos < str.length(); ++pos) {    // collect the bases from the string
        if (str[pos] != 'A' && str[pos] != 'C' && str[pos] != 'G' && str[pos] != 'T') {
            str = str.substr(pos+1, string::npos);
            goto next_kmer;    // unknown base, start a new k-mer from the beginning
        }
        kmer::shift_right(kmer, str[pos]);    // shift each base into the bit sequence

        if (pos+1 >= kmer::k) {
            rcmer = kmer;
            if (reverse) kmer::reverse_complement(rcmer, true);    // invert the k-mer, if necessary
            color::set(kmer_table[rcmer], color);    // update the k-mer with the current color
        }
    }
}

/**
 * This function extracts k-mers from a sequence and adds them to the hash table.
 *
 * @param str dna sequence
 * @param color color flag
 * @param reverse merge complements
 * @param max_iupac allowed number of ambiguous k-mers per position
 */
void graph::add_kmers(string& str, uint64_t& color, bool& reverse, uint64_t& max_iupac) {

    unordered_set<kmer_t> ping;    // create a new empty set for the k-mers
    unordered_set<kmer_t> pong;    // create another new set for the k-mers
    bool ball; bool wait;    // indicates which of the two sets should be used

    vector<uint8_t> factors;    // stores the multiplicity of iupac bases
    long double product;    // stores the overall multiplicity of the k-mers
    uint64_t pos;    // current position in the string, from 0 to length

    kmer_t kmer0;    // create an empty bit sequence for the initial k-mer
    kmer_t rcmer;    // create a bit sequence for the reverse complement

next_kmer:
    ping.clear(); pong.clear(); factors.clear();
    if (str.length() < kmer::k) {
        return;    // not enough characters for a k-mer, abort
    }
    ball = true; wait = false; product = 1; pos = 0;
    (ball ? ping : pong).emplace(kmer0);

    for (; pos < str.length(); ++pos) {    // collect the bases from the string
        if (str[pos] == '.' || str[pos] == '-') {
            str = str.substr(pos+1, string::npos);
            goto next_kmer;    // gap character, start a new k-mer from the beginning
        }
        iupac_calc(product, factors, str[pos]);

        if (product <= max_iupac) {    // check if there are too many ambiguous k-mers
            if (wait) {
                str = str.substr(pos-kmer::k+1, string::npos);
                goto next_kmer;    // start a new k-mer from the beginning
            }
            iupac_shift(ball ? ping : pong, !ball ? ping : pong, str[pos]);
            ball = !ball;    // shift each base in, resolve iupac character
        } else { wait = true; continue; }

        if (pos+1 >= kmer::k) {
            for (auto& kmer : (ball ? ping : pong)) {    // iterate over the current set of ambiguous k-mers
                rcmer = kmer;
                if (reverse) kmer::reverse_complement(rcmer, true);    // invert the k-mer, if necessary
                color::set(kmer_table[rcmer], color);    // update the k-mer with the current color
            }
        }
    }
}

/**
 * This function calculates the multiplicity of iupac k-mers.
 *
 * @param product overall multiplicity
 * @param factors per base multiplicity
 * @param input iupac character
 */
void graph::iupac_calc(long double& product, vector<uint8_t>& factors, char& input) {

    switch (input) {
        case 'A': case 'C': case 'G': case 'T':
            product *= 1;
            factors.emplace_back(1);
    }
    switch (input) {
        case 'R': case 'Y': case 'S': case 'W': case 'K': case 'M':
            product *= 2;
            factors.emplace_back(2);
    }
    switch (input) {
        case 'B': case 'D': case 'H': case 'V':
            product *= 3;
            factors.emplace_back(3);
    }
    switch (input) {
        case 'N':
            product *= 4;
            factors.emplace_back(4);
    }
    if (factors.size() > kmer::k) {
        product /= *factors.begin();
        factors.erase(factors.begin());
    }
}

/**
 * This function shifts a base into a set of ambiguous iupac k-mers.
 *
 * @param prev set of k-mers
 * @param next set of k-mers
 * @param input iupac character
 */
void graph::iupac_shift(unordered_set<kmer_t>& prev, unordered_set<kmer_t>& next, char& input) {

    kmer_t temp; char base;
    while (!prev.empty()) {    // extend each previous k-mer
        switch (input) {
            case 'A': case 'R': case 'W': case 'M':
            case 'D': case 'H': case 'V': case 'N':
                temp = *prev.begin(); base = 'A';
                kmer::shift_right(temp, base);
                next.emplace(temp);
        }
        switch (input) {
            case 'C': case 'Y': case 'S': case 'M':
            case 'B': case 'H': case 'V': case 'N':
                temp = *prev.begin(); base = 'C';
                kmer::shift_right(temp, base);
                next.emplace(temp);
        }
        switch (input) {
            case 'G': case 'R': case 'S': case 'K':
            case 'B': case 'D': case 'V': case 'N':
                temp = *prev.begin(); base = 'G';
                kmer::shift_right(temp, base);
                next.emplace(temp);
        }
        switch (input) {
            case 'T': case 'Y': case 'W': case 'K':
            case 'B': case 'D': case 'H': case 'N':
                temp = *prev.begin(); base = 'T';
                kmer::shift_right(temp, base);
                next.emplace(temp);
        }
        prev.erase(prev.begin());
    }
}

/**
 * This function iterates over the hash table and calculates the split weights.
 *
 * @param mean weight function
 * @param verbose print progress
 */
void graph::add_weights(double mean(uint32_t&, uint32_t&), bool& verbose) {

    double min_value = numeric_limits<double>::min();    // current min. weight in the top list (>0)
    uint64_t cur = 0, prog = 0, next;
    uint64_t max = kmer_table.size();
loop:
    for (auto& it : kmer_table) {    // iterate over k-mer hash table
        if (verbose) {
            next = 100*cur/max;
             if (prog < next)  cout << "\33[2K\r" << "Processing splits... " << next << "%" << flush;
            prog = next; cur++;
        }
        color_t& color = it.second;    // get the color set for each k-mer
        bool pos = color::complement(color, true);    // invert the color set, if necessary
        if (color == 0) continue;    // ignore empty splits
        array<uint32_t,2>& weight = color_table[color];    // get the weight and inverse weight for the color set

        double old_value = mean(weight[0], weight[1]);    // calculate the old mean value
        if (old_value >= min_value) {    // if it is greater than the min. value, find it in the top list
            auto range = split_list.equal_range(old_value);    // get all color sets with the given weight
            for (auto it = range.first; it != range.second; ++it) {
                if (it->second == color) {    // iterate over the color sets to find the correct one
                    split_list.erase(it);    // erase the entry with the old weight
                    break;
                }
            }
        }
        weight[pos]++;    // update the weight or the inverse weight of the current color set

        double new_value = mean(weight[0], weight[1]);    // calculate the new mean value
        if (new_value >= min_value) {    // if it is greater than the min. value, add it to the top list
            split_list.emplace(new_value, color);    // insert it at the correct position ordered by weight
            if (split_list.size() > t) {
                split_list.erase(--split_list.end());    // if the top list exceeds its limit, erase the last entry
                min_value = split_list.rbegin()->first;    // update the min. value for the next iteration
            }
        }
    }
}

/**
 * This function adds a single split (weight and colors) to the output list.
 *
 * @param weight split weight
 * @param color split colors
 */
void graph::add_split(double& weight, color_t& color) {
    split_list.emplace(weight, color);    // insert it at the correct position ordered by weight
    if (split_list.size() > t) {
        split_list.erase(--split_list.end());    // if the top list exceeds its limit, erase the last entry
    }
}

/**
 * This function filters a greedy maximum weight tree compatible subset.
 *
 * @param verbose print progress
 */
void graph::filter_strict(bool& verbose) {
	filter_strict(0, verbose);
}


/**
 * This function filters a greedy maximum weight tree compatible subset and returns a newick string.
 *
 * @param map function that maps an integer to the original id, or null if no newick output wanted.
 * @param verbose print progress
 */
string graph::filter_strict(std::function<string(const uint64_t&)> map, bool& verbose) {
    auto tree = vector<color_t>();    // create a set for compatible splits
    auto it = split_list.begin();
    uint64_t cur = 0, prog = 0, next;
    uint64_t max = split_list.size();
loop:
    while (it != split_list.end()) {
        if (verbose) {
            next = 100*cur/max;
             if (prog < next)  cout << "\33[2K\r" << "Filtering splits... " << next << "%" << flush;
            prog = next; cur++;
        }
        if (test_strict(it->second, tree)) {
            tree.emplace_back(it->second);
            ++it; goto loop;    // if compatible, add the new split to the set
        }
        it = split_list.erase(it);    // otherwise, remove split
    }
    if(map){
		set * root = build_tree(tree);
		return print_tree(root,map) + ";\n";
	}else{
		return "";
	}
}



/**
 * This function filters a greedy maximum weight weakly compatible subset.
 *
 * @param verbose print progress
 */
void graph::filter_weakly(bool& verbose) {
    auto network = vector<color_t>();    // create a set for compatible splits
    auto it = split_list.begin();
    uint64_t cur = 0, prog = 0, next;
    uint64_t max = split_list.size();
loop:
    while (it != split_list.end()) {
        if (verbose) {
            next = 100*(cur*cur)/(max*max);
             if (prog < next)  cout << "\33[2K\r" << "Filtering splits... " << next << "%" << flush;
            prog = next; cur++;
        }
        if (test_weakly(it->second, network)) {
            network.emplace_back(it->second);
            ++it; goto loop;    // if compatible, add the new split to the set
        }
        it = split_list.erase(it);    // otherwise, remove split
    }
}




/**
 * This function filters a greedy maximum weight n-tree compatible subset.
 *
 * @param n number of trees
 * @param verbose print progress
 */
void graph::filter_n_tree(uint64_t n, bool& verbose) {
    filter_n_tree(n, 0,  verbose);
}

/**
 * This function filters a greedy maximum weight n-tree compatible subset and returns a string with all trees in newick format.
 *
 * @param n number of trees
 * @param map function that maps an integer to the original id, or null
 * @param verbose print progress
 */
string graph::filter_n_tree(uint64_t n, std::function<string(const uint64_t&)> map, bool& verbose) {
    auto forest = vector<vector<color_t>>(n);    // create a set for compatible splits
    auto it = split_list.begin();
    uint64_t cur = 0, prog = 0, next;
    uint64_t max = split_list.size();
loop:
    while (it != split_list.end()) {
        if (verbose) {
            next = 100*cur/max;
             if (prog < next)  cout << "\33[2K\r" << "Filtering splits... " << next << "%" << flush;
            prog = next; cur++;
        }
       for (auto& tree : forest)
        if (test_strict(it->second, tree)) {
            tree.emplace_back(it->second);
            ++it; goto loop;    // if compatible, add the new split to the set
        }
        it = split_list.erase(it);    // otherwise, remove split
    }
    // output
    string s="";
	if(map){
		for (auto& tree : forest) {
			set * root = build_tree(tree);
			s+=print_tree(root,map) + ";\n";
		}
	}
	return s;
		
}





/**
 * This function tests if a split is compatible with an existing set of splits.
 *
 * @param color new split
 * @param color_set set of splits
 * @return true, if compatible
 */
bool graph::test_strict(color_t& color, vector<color_t>& color_set) {
    for (auto& elem : color_set) {
        if (!color::is_compatible(elem, color)) {
            return false;    // compare to each split in the set
        }
    }
    return true;
}

/**
 * This function tests if a split is weakly compatible with an existing set of splits.
 *
 * @param color new split
 * @param color_set set of splits
 * @return true, if weakly compatible
 */
bool graph::test_weakly(color_t& color, vector<color_t>& color_set) {
    for (auto& elem1 : color_set) {
        for (auto& elem2 : color_set) {
            if (elem1 != elem2) {
                if (!color::is_weakly_compatible(elem1, elem2, color)) {
                    return false;    // compare to each split in the set
                }
            }
        }
    }
    return true;
}


/**
 * This function recursively refines a given set/tree structure by a given split
 *
 * @param current_set node of currently considered (sub-)set/tree structure
 * @param split color set to refine by
 * @return whether or not the given split is compatible with the set/tree structure
 */
bool graph::refine_tree(set * current_set, color_t& split, color_t& allTaxa){
     
    // possible cases:
    // splitsize <2: nothing has to be done
    // split equals one subset -> warning: split twice
    // split is fully contained in one subset -> recurse
    // inverse split ... (i.e. split covers one subset partially) -> recurse with inverse
    // split covers several subsets completely -> introduce new split
    
	 vector<set *> * subsets= &current_set->subsets;
	 color_t taxa=current_set->taxa;
	 
    if (color::size(split,false)<2 or color::size(allTaxa,false)-color::size(split,false)<2){ return true; }

	vector<set *> fullycoveredsubsets={};
    set * partiallycoveredsubset = nullptr;

 	for (set * subset : *subsets) {
		color_t subtaxa = subset->taxa;
		if (split == subtaxa){
			return true;
		}
		//split.issubset(subtaxa)?
		if( (split & subtaxa) == split){ return refine_tree(subset,split,allTaxa); }
		// subtaxa.issubset(split):
		if( (subtaxa & split) == subtaxa){ fullycoveredsubsets.push_back(subset); }
		// elif not subtaxa.isdisjoint(split): # does intersect
		else if((subtaxa & split) != 0b0u){
			// if partiallycoveredsubset:
			if(partiallycoveredsubset!=nullptr){ return false; }//there cannot be more than one
			else{ partiallycoveredsubset=subset; }
		}
 	}            
    if (partiallycoveredsubset!=nullptr){
		if(fullycoveredsubsets.size()==subsets->size()-1){
			//recurse into this subset with inverse split
			color_t inversesplit=color::complement(split,false);
			//  if inversesplit.issubset(partiallycoveredsubset[1]):
			if ( (inversesplit & partiallycoveredsubset->taxa ) == inversesplit ){
				return refine_tree(partiallycoveredsubset,inversesplit,allTaxa);
			} else { return false; }
		} else {  return false; }
	} else if (fullycoveredsubsets.size()>1){
		// introduce new split
        color_t newsubtaxa=0b0u;
		for(set * subset : fullycoveredsubsets){ newsubtaxa|=subset->taxa; }
		//get weight of split
		double weight = 0;
		auto it = split_list.begin();
		while (it != split_list.end()) {
			if (it->second==split){
				weight=it->first;
				break;
			}
			it++;
		}		
		set * newset = newSet(newsubtaxa,weight,fullycoveredsubsets);
		// remove old sets
        for(set * subset : fullycoveredsubsets){
            //subsets.remove(subset)
 			subsets->erase(std::remove(subsets->begin(), subsets->end(), subset), subsets->end());
		}
		// add new set
		subsets->push_back(newset);
		return true;
	} else {
        std::cerr << "ERROR: this cannot be: just one fully covered subset and nothing else!?" << "\n" << flush;
        exit(EXIT_FAILURE);
	}
	std::cerr << "ERROR: This point in code should never be reached" << "\n" << flush;
    exit(EXIT_FAILURE);
	
 }
 
 

/**
 * This function returns a tree structure (struct set) generated from the given list of color sets
 *
 * @param color_set list of color sets
 * @return tree structure (struct set)
 */
set * graph::build_tree(vector<color_t>& color_set){

	uint64_t n = color::n;

	//initialize set of trivial splits
	vector<set *> subsets = {};
	color_t allTaxa=0b0u;
	for (uint64_t i = 0; i < n; i++) {
		color_t leaf=0b0u;
		color::set(leaf,i);
		color::set(allTaxa,i);
		vector<set *> emptyset = {};
		// get weight
		double weight = 0;
		auto it = split_list.begin();
		while (it != split_list.end()) {
			if (it->second==leaf){
				weight=it->first;
				break;
			}
			it++;
		}
		set * newset = newSet(leaf,weight,emptyset);
		subsets.push_back(newset);
	}
	set * sets = newSet(allTaxa,0,subsets);

	
	for (color_t split : color_set){
		//split if possible
		if (!refine_tree(sets,split,allTaxa)){
			std::cerr << "ERROR: splits are incompatible\n" << flush;
			exit(EXIT_FAILURE);
		}
	}

	return sets;
}



/**
* This function returns a newick string generated from the given tree structure (set)
*
* @param root root of the tree/set structure
* @return newick string
*/
string graph::print_tree(set * root, std::function<string(const uint64_t&)> map){
	 vector<set *> subsets=root->subsets;
	 color_t taxa=root->taxa;
    if (subsets.size()==0){ // leaf set
        if (color::size(taxa,false)==0) {
            std::cerr << "ERROR: child with no taxon!?\n" << flush;
			exit(EXIT_FAILURE);
		} else if(color::size(taxa,false)==1){
                return map(color::pos(taxa))+":"+to_string(root->weight);
		} else {
            std::cerr << "ERROR: child with more than one taxon!?\n" << flush;
			exit(EXIT_FAILURE);
		}
	}
    else {
        string s="(";
		for (set * subset : subsets){
			s += print_tree(subset,map);
			if (subset!=subsets.back()) { s += ","; }
		}
        s += "):";
		s+=to_string(root->weight);
        return s;
	}
}





 
 
 
 
 
 
 
 
 
 
