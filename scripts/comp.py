#!/usr/bin/env python3

from __future__ import print_function
import sys
import os
import pygtrie


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def unify(split, taxa):
    ls=len(split)
    lt=len(taxa)
    if ls==(lt-ls):
        sorted_inv = sorted([item for item in taxa if item not in split])
        sorted_split=sorted(split)
        if sorted_inv<sorted_split:
            return sorted_inv
        else:
            return sorted_split
    if (lt-ls)<ls:
        inv = [item for item in taxa if item not in split]
        return sorted(inv)
    else:
        return sorted(split)

def readtaxa(filename):
    taxa=set()
    for line in (s.strip() for s in open(filename)):
        taxa.add(line)
    return(taxa)


def readfile(filename,taxa,add_weight):
    splits=pygtrie.StringTrie()
    for line in (s.strip() for s in open(filename)):
        #split line
        fields = line.split('\t')
        weight = float(fields[0])+add_weight
        split = fields[1:]
        #store tax names if necessary
        for f in split:
            if f not in taxa:
                eprint("taxa from split not in taxa: "+f)
                exit(1)
        split=unify(split,taxa)
        split_txt="/".join(split)
        #update or store split
        if split_txt in splits:
            splits[split_txt]+=weight
        else:
            splits[split_txt]=weight

    return splits




def get_top(splitlist,l,num): #l=len(splitlist)
    

    # too short to be filtered at all?
    if l<=num:
        return splitlist
    # else..
    
    #sort by weights
    both_sorted=sorted([(w,s) for (s,w) in splitlist.items()], reverse=True)
    
    #cut top num
    for (_,s) in both_sorted[num:]:
        del splitlist[s]



def remove_leaves(splitlist,taxa):
    for t in taxa:
        if t in splitlist:
            del splitlist[t]


def split_comp(donor,reference):
    
    #sort by weights
    donor_sorted=sorted([(w,s) for (s,w) in donor.items()], reverse=True)
    
    ref_list=reference.keys()
    n=len(ref_list)
    
    w_all=0.
    w_corr=0.
    num_all=0
    num_corr=0
    for (w,s) in donor_sorted:
        w_all+=w
        num_all+=1
        if s in ref_list:
            ref_list.remove(s)
            w_corr+=w
            num_corr+=1
            #print("\t".join([str(num_all),str(w),"-1",str(len(s.split('/')))]))
        #else:
            #print("\t".join([str(num_all),"-1",str(w),str(len(s.split('/')))]))
#            eprint(s.split('/'))
        #if len(ref_list) == 0:
            #break
    #for r in ref_list:
        #print("\t".join(["-1","-1","-1",str(len(r.split('/')))]))
    #eprint(w_all,w_corr,num_all,num_corr,w)

    print("#precision\trecall")
    precision=(1.0*num_corr)/(1.0*num_all)
    recall=(1.0*num_corr)/(1.0*n)
    print("\t".join([str(precision),str(recall)]))
                     
    return w_corr

                        


    

if len(sys.argv)>4:
    add_weight=float(sys.argv[4])
else:
    add_weight=0


if len(sys.argv)<4:
    eprint("Usage: sans2nexus.py <list of splits 1> <list of splits 2> <list of genomes (file names)>")
    eprint("list of splits in SANS output format: Per split one line: weight genomeA genomeB ...")
    eprint("Output on stdout.")
    sys.exit(1)


   
eprint("read taxa")
taxa=readtaxa(sys.argv[3])
eprint(len(taxa))

eprint("read split file 1")
splitlist1=readfile(sys.argv[1],taxa,add_weight)
eprint("read split file 2")
splitlist2=readfile(sys.argv[2],taxa,add_weight)
eprint("found:")
l1=len(splitlist1)
l2=len(splitlist2)
eprint(l1,l2)

#eprint("get copy without leaves")
#leaflesssplits1=splitlist1.copy()
#remove_leaves(leaflesssplits1,taxa)
#leaflesssplits2=splitlist2.copy()
#remove_leaves(leaflesssplits2,taxa)
#eprint("left:")
#ll1=len(leaflesssplits1)
#ll2=len(leaflesssplits2)
#eprint(ll1,ll2)

#eprint("filter top n**2-n ("+str(len(taxa)**2-len(taxa))+") from leafless sets")
#get_top(leaflesssplits1,ll1,len(taxa)**2-len(taxa))
#get_top(leaflesssplits2,ll2,len(taxa)**2-len(taxa))

#eprint("compute RF on leafless sets")
#rf_l=comp_diff_sqr(leaflesssplits1,leaflesssplits2,splits2lists(leaflesssplits1,leaflesssplits2))
#eprint("compute FM on leafless sets")
#fm_l=comp_diff_sqr(leaflesssplits1,leaflesssplits2,get_pairwise_ds(leaflesssplits1,leaflesssplits2,taxa))


#eprint("filter top n**2 ("+str(len(taxa)**2)+") from complete sets")
#get_top(splitlist1,l1,len(taxa)**2)
#get_top(splitlist2,l2,len(taxa)**2)

#eprint("compute RF")
#rf=comp_diff_sqr(splitlist1,splitlist2,splits2lists(splitlist1,splitlist2))
#eprint("compute FM")
#fm=comp_diff_sqr(splitlist1,splitlist2,get_pairwise_ds(splitlist1,splitlist2,taxa))
eprint("compute precision and recall:")
sc=split_comp(splitlist1,splitlist2)


# output
#print("\t".join(map(str,[rf,rf_l,fm,fm_l,sc])))
#print("\t".join([str(diff_taxa)]))


