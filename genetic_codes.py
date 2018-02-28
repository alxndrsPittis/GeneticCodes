import cPickle
from pandas import DataFrame
from collections import Counter, defaultdict
import os

from ete3 import SeqGroup

if os.path.exists('alg.pkl'):
    print 'loading from pkl file'
    a = cPickle.load(open('alg.pkl'))
else: 
    alg = SeqGroup('formatted_MG_seqs.faa.final_tree.fa.gz')
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
for col in a:
    counter = Counter(a[col])
    refaa, num = counter.most_common(1)[0]
    frac = num/nsp
    if refaa != '-':
        print col, valid_cols, refaa, frac
    
    if refaa != '-' and frac > 0.75:
        valid_cols += 1
        variants = a[a[col] != refaa][col].to_dict()
        for sp, var in variants.iteritems():
            spvariants[sp].update([(refaa, var)])

    if valid_cols > 500:
        break

            
for sp, varcounter in spvariants.iteritems():
    if varcounter.most_common(1)[0][1] > 10:
        print sp, ncbi.translate_to_names([int(sp.split(".")[0])])[0], varcounter.most_common(5)


        

    
# for col in a:
#     print a[col]
#     raw_input()

