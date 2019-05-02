#!/usr/bin/env python3

import dendropy
import sys

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



if len(sys.argv)<2:
    eprint("Usage: newick2sans.py <newick file> [<list of genomes (file names) to be considered, i.e., prune tree to given list>]")
    eprint("Output in SANS output format on stdout.")
    sys.exit(1)


tree = dendropy.Tree.get(
        path=sys.argv[1],
        schema="newick",
        preserve_underscores=True)
eprint(len(tree.leaf_nodes()))


#if len(sys.argv)>2:
    ## read mapping
    #mapping=dict()
    #inv_mapping=dict()
    #for line in (s.strip() for s in open(sys.argv[2])):
        #items=line.split('\t')
        #mapping[items[0]]=items[1]
        #inv_mapping[items[1]]=items[0]
    ##prune
    #node_filter_fn = lambda nd: nd.is_internal() or nd.taxon.label in mapping.keys()
    #tree = tree.extract_tree(node_filter_fn=node_filter_fn)
    #eprint(len(tree.leaf_nodes()))

if len(sys.argv)>2:
    eprint("pruning tree...")
    leaflables=[nd.taxon.label for nd in tree.leaf_nodes()]
    #read taxa
    taxa=set()
    for line in (s.strip() for s in open(sys.argv[2])):
        taxa.add(line)
        #if line not in mapping.values():
            #eprint("taxa from list not in mapping (tree pruned) : "+line)
        #elif inv_mapping[line] not in leaflables:
            #eprint("taxa from list not in tree (tree pruned) : "+line+" <- "+inv_mapping[line])
        if line not in leaflables:
            eprint("taxa from list not in tree: "+line)
    eprint("taxa read: "+str(len(taxa)))
#    for l in leaflables:
#        if  l not in taxa:
#            eprint("taxa from tree not in taxa list (tree pruned) : "+l)
    node_filter_fn = lambda nd: nd.is_internal() or nd.taxon.label in taxa
    tree = tree.extract_tree(node_filter_fn=node_filter_fn)
    #eprint(len(tree.leaf_nodes()))


tree.encode_bipartitions(suppress_unifurcations=True,collapse_unrooted_basal_bifurcation=True)
for edge in tree.preorder_edge_iter():
#for node in tree:
    if edge.length:
        split=str(edge.length)+"\t"
    else:
        split="0.0\t"
#    if len(sys.argv)==2:
    split+="\t".join([t.label for t in edge.bipartition.leafset_taxa(tree.taxon_namespace)])
#    else:
#        split+="\t".join(map(map_taxa,[t.label for t in edge.bipartition.leafset_taxa(tree.taxon_namespace)]))
    print(split)
