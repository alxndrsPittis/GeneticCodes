import sys
import cPickle
from pandas import DataFrame
from collections import Counter, defaultdict
import os
from glob import glob
from sys import exit
from ete3 import SeqGroup
from ete3 import NCBITaxa

ncbi = NCBITaxa()

path = sys.argv[1] + "*clustalo"
infiles = glob(path)
print "%s infiles" % len(infiles)

valid_cols = 0
spvariants = defaultdict(Counter)
refaas = []

for infile in infiles:
#for infile in glob("formatted_MG_seqs.faa.final_tree.fa"):
    print infile
    if os.stat(infile).st_size == 0:
        continue
    alg = SeqGroup(infile)
    alg_matrix = []
    labels = []    
    
    for name, seq, _ in alg:
        # Replace trailing gaps with # and * for stop
        # Count end positions
        for n, aa in enumerate(seq[::-1]):
            if aa != "-":
                break
            
        # Seems that is not possible to know where the stop codon would align..
        # seq = seq[:len(seq)-n] + "*" + "#" * (n-1)
        # So just keep track of trailing gaps, as compared to internal
        seq = seq[:len(seq)-n] + "*" * n
                        
        aas = [aa for aa in seq.upper()]
        alg_matrix.append(aas)
        labels.append(name)

    a = DataFrame(alg_matrix, labels)
    nsp = float(len(a))

    for col in a:
        counter = Counter(a[col])
        refaa, num = counter.most_common(1)[0]
        frac = num/nsp

        #if refaa != '-':
            #pass
            #print col, valid_cols, refaa, frac

        if refaa !=  '-' and frac > 0.66:
            valid_cols += 1
            variants = a[a[col] != refaa][col].to_dict()
            refaas.append(refaa)
            for sp, var in variants.iteritems():
                sp = sp.split(".")[0]
                spvariants[sp].update([(refaa, var)])
    
        #if valid_cols > 500:
        #    break

refaacounter = Counter(refaas)
            
for sp, varcounter in spvariants.iteritems():
    try:
        sp_name = ncbi.translate_to_names([int(sp.split(".")[0])])[0]
    except ValueError:
        sp_name = "oxymonad-%s" % sp
    for varc in varcounter:
        ratio = varcounter[varc] / float(refaacounter[varc[0]])
        if ratio > 0.25:
            print sp, sp_name, "%s\t%s/%s" % ( "->".join(varc), varcounter[varc], refaacounter[varc[0]] )
    # #print varcounter.most_common(1)[0], varcounter[('W', 'X')]
    # most_common = varcounter.most_common(1)[0]
    # ratio = (most_common[1] / float(refaacounter[most_common[0][0]]))
    # #if varcounter.most_common(1)[0][1] > 10 and "-" not in varcounter.most_common(1)[0][0]:
    # if ratio > 0.33:# and "-" not in most_common[0]:
    #     print sp, ncbi.translate_to_names([int(sp.split(".")[0])])[0]
    #     for varc in varcounter:
    #         #if varcounter[varc] > 2:#/float(refaacounter[varc[0]]) > 0.2:
    #         try:
    #             print "%s\t%s/%s" % ( "->".join(varc), varcounter[varc], refaacounter[varc[0]] )
    #         except TypeError:
    #             print "%s\t%s/%s" % ( "%s->None" % varc[0], varcounter[varc], refaacounter[varc[0]] )
    #     #print sp, varcounter.most_common(5)
        

    
# for col in a:
#     print a[col]
#     raw_input()

