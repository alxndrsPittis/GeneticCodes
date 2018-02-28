import cPickle
from pandas import DataFrame
from collections import Counter, defaultdict
import os

from ete3 import SeqGroup

if os.path.exists('alg.pkl'):
    print 'loading from pkl file'
    a = cPickle.load(open('alg.pkl'))
else: 
    alg = SeqGroup('formatted_MG_seqs.faa.final_tree.fa')
    #alg = SeqGroup('Burki_first10.aln.fa')
    #alg = SeqGroup('test_52.fa')
    alg_matrix = []
    labels = []
    for name, seq, _ in alg:
        aas = [aa for aa in seq.upper()]
        alg_matrix.append(aas)
        labels.append(name)

    a = DataFrame(alg_matrix, labels)
    cPickle.dump(a, open('alg.pkl', 'wb'), 2)

nsp = float(len(a))
spvariants = defaultdict(Counter)
valid_cols = 0
refaas = []
for col in a:
    counter = Counter(a[col])
    refaa, num = counter.most_common(1)[0]
    frac = num/nsp
    if refaa != '-':
        print col, valid_cols, refaa, frac
    
    if refaa != '-' and frac > 0.5:
        valid_cols += 1
        variants = a[a[col] != refaa][col].to_dict()
        refaas.append(refaa)
        for sp, var in variants.iteritems():
            spvariants[sp].update([(refaa, var)])

    if valid_cols > 500:
        break

refaacounter = Counter(refaas)

for sp, varcounter in spvariants.iteritems():
    #print varcounter.most_common(1)[0], varcounter[('W', 'X')]
    most_common = varcounter.most_common(1)[0]
    ratio = (most_common[1] / float(refaacounter[most_common[0][0]]))
    #if varcounter.most_common(1)[0][1] > 10 and "-" not in varcounter.most_common(1)[0][0]:
    if ratio > 0.33 and "-" not in most_common[0]:    
        print sp
        for varc in varcounter:
            if varcounter[varc]/float(refaacounter[varc[0]]) > 0.2:
                try:
                    print "%s\t%s/%s" % ( "->".join(varc), varcounter[varc], refaacounter[varc[0]] )
                except TypeError:
                    print "%s\t%s/%s" % ( "%s->None" % varc[0], varcounter[varc], refaacounter[varc[0]] )
        #print sp, varcounter.most_common(5)


        

    
# for col in a:
#     print a[col]
#     raw_input()

