#!/usr/bin/env python3


from __future__ import print_function
import sys
import os
import pygtrie
from math import log
from math import sqrt

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



def readfile(filename,taxa):
    splits=pygtrie.StringTrie()
    for line in (s.strip() for s in open(filename)):
        #split line
        fields = line.split('\t')
        weight = float(fields[0])
        split_orig = fields[1:]
        #store tax names if necessary
        split=[] #storing ids
        for f in split_orig:
            if f not in taxa.keys():
                eprint("taxa from split not in taxa: "+f)
                exit(1)
            else:
                split.append(taxa[f])
        if len(split)==len(taxa) or len(split)==0:
            continue
        split=unify(split,taxa.values())
        split_txt="/".join(map(str,split))
        #update or store split
        if split_txt in splits:
            splits[split_txt]+=weight
        else:
            splits[split_txt]=weight

    return splits


                        

if len(sys.argv)>3:
    min_weight=int(sys.argv[3])
else:
    min_weight=-1

if len(sys.argv)>4:
    add_weight=int(sys.argv[4])
else:
    add_weight=0



def readtaxa(filename):
    taxa=dict()
    i=1
    for line in (s.strip() for s in open(filename)):
        taxa[line]=i
        i+=1
    return(taxa)


if len(sys.argv)<3:
    eprint("Usage: sans2nexus.py <sans output file> <list of genomes (file names)>")
    eprint("Output in NEXUS format on stdout.")
    sys.exit(1)


eprint("read taxa")
taxa=readtaxa(sys.argv[2])
ntax=len(taxa)
eprint(ntax)

eprint("read split file")
splitlist=readfile(sys.argv[1],taxa)
eprint("read:")
l=len(splitlist)
eprint(l)

#eprint("filter top 10*n, i.e. "+str(ntax*10))
#sorted_splits=sorted([(w,s) for (s,w) in splitlist.iteritems()], reverse=True)
#splitlist=sorted_splits[:(ntax*10)]

   
# output nexus format stuff
print("#nexus\n\nBEGIN Taxa;\nDIMENSIONS ntax=%s;\nTAXLABELS"%(ntax))

# output mapping filename->taxid
for tupel in sorted( ((v,k) for k,v in taxa.items())):
    print("[%s] '%s' "%(tupel))

# prepare output of splits and count
nsplits=0
splitstring=""
for (weight,split_str) in splitlist:
    if split_str=="":
        split=[]
    else:
        split=split_str.split('/')
    if (weight>min_weight or len(split)==1 or len(split)==ntax-1) and len(split)>0: #don't output split all versus none and too small splits
        nsplits+=1
        #splitstring+="%s\t%s\n"%(len(split),weight)
        #splitstring+="[%s, size=%s]\t%s\t%s,\n"%(nsplits,len(split),log(weight)*len(split),' '.join(map(str,split)))
        splitstring+="[%s, size=%s]\t%s\t%s,\n"%(nsplits,len(split),weight,' '.join(map(str,split)))
        #splitstring+="[%s, size=%s]\t%s\t%s,\n"%(nsplits,len(split),1.0*weight/sqrt(len(split)*(ntax-len(split))),' '.join(map(str,split)))
        #splitstring+="[%s, size=%s]\t%s\t%s,\n"%(nsplits,len(split),1.0*weight/len(split),' '.join(map(str,split)))
        #splitstring+="[%s, size=%s]\t%s\t%s,\n"%(nsplits,len(split),weight/100.0,' '.join(map(str,split)))
        #splitstring+="[%s, size=%s]\t%s\t%s,\n"%(nsplits,len(split),log(weight),' '.join(map(str,split)))


# output nexus format stuff
print(";\nEND; [Taxa]\n\n\n\nBEGIN Splits;\nDIMENSIONS ntax=%s nsplits=%s;\nMATRIX"%(ntax,nsplits))

#output splits
print(splitstring)

# output nexus format stuff
print(";\nEND; [Splits]")


print("BEGIN st_Assumptions;\nuptodate;\nsplitstransform=EqualAngle;\nSplitsPostProcess filter=none;\nexclude  no missing;\nautolayoutnodelabels;\nEND; [st_Assumptions]")

eprint("output:")
eprint(nsplits)