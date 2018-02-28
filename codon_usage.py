import sys
import os
import cPickle
from pandas import DataFrame
from collections import Counter, defaultdict
import itertools
from sys import exit
from ete3 import SeqGroup
from ete3 import NCBITaxa
from Bio.Data import CodonTable

ambiguous_map = { "Y" : ["C", "T"],
                  "R" : ["A", "G"],
                  "W" : ["A", "T"],
                  "S" : ["G", "C"],
                  "K" : ["T", "G"],
                  "M" : ["C", "A"],
                  "D" : ["A", "G", "T"],
                  "V" : ["A", "C", "G"],
                  "H" : ["A", "C", "T"],
                  "B" : ["C", "G", "T"],
                  "N" : ["A", "G", "C", "T"]
    }

def ambiguity_base(base_counter):
    # distribute ambiguous base frequencies
    base_contribution = defaultdict(dict)
    ambiguous = set(base_counter.keys()) - set(bases)
    if ambiguous:
        for amb in ambiguous:
            # contribution of Y to C : Cn/(Cn+Tn)
            nreps = float(sum([base_counter[base] for base in ambiguous_map[amb]]))
            for base in ambiguous_map[amb]:
                base_contribution[amb][base] = round(base_counter[base] / nreps, 2)
        # add contributions and remove ambiguous base
        for amb in base_contribution:
            for base in base_contribution[amb]:
                # total contribution of Y to C : Yn*Cn/(Cn+Tn)
                base_counter[base] += round(base_counter[amb] * base_contribution[amb][base], 2)
                del base_counter[amb]
        return base_contribution
    else:
        return False
            
def ambiguity_codon(codon_counter, base_contribution):
    # distribute ambiguous codons frequencies based on bases
    for triplet in codon_counter.keys():
        # build iterables for product and keep indices
        combinations = []
        amb_indices = []
        for i, base in enumerate(triplet):
            if base in ambiguous_map:
                combinations.append(ambiguous_map[base])
                amb_indices.append(i)
            else:
                combinations.append(base)
        if amb_indices:
            for combination in itertools.product(*combinations):
                codon = "".join(combination)
                # multiply single base contributions, based on indices kept
                contribution = reduce(lambda x, y: x*y, [base_contribution[triplet[i]][codon[i]] for i in amb_indices])
                # total contribution of codon is a product of its count
                codon_counter[codon] += codon_counter[triplet]*contribution
            del codon_counter[triplet]

def translate(seq):
    seq = seq.lower().replace('\n', '').replace(' ', '')
    peptide = ''

    for i in xrange(0, len(seq), 3):
        codon = seq[i: i+3]
        amino_acid = codon_table.get(codon, '*')
        if amino_acid != '*':
            peptide += amino_acid
        else:
            break
        
    return peptide

infile = sys.argv[1]
outpath = sys.argv[2]


# store all possible start and stop codons
all_stops = []
all_starts = []
for code_id in CodonTable.unambiguous_dna_by_id:
    all_stops += CodonTable.unambiguous_dna_by_id[code_id].stop_codons
    all_starts += CodonTable.unambiguous_dna_by_id[code_id].start_codons
all_stops = set(all_stops)
all_starts = set(all_starts)

# get all codons
uni_table = CodonTable.unambiguous_dna_by_id[1]
all_codons = (uni_table.stop_codons + uni_table.forward_table.keys())
bases = ['T', 'C', 'A', 'G']
# sort according to table
all_codons.sort(key=lambda x: (bases.index(x[0]), bases.index(x[1]), bases.index(x[2])))

if os.path.exists(outpath + '.spcodons.pkl'):
    print 'loading counts from pkl file'
    cdna = cPickle.load(open(outpath + '.alg.pkl'))
    spcodons = cPickle.load(open(outpath + '.spcodons.pkl'))
    spbases = cPickle.load(open(outpath + '.spbases.pkl'))
    spstarts = cPickle.load(open(outpath + '.spstarts.pkl'))
    spstops = cPickle.load(open(outpath + '.spstops.pkl'))
    print 'loading counts done'
else:
    if os.path.exists(outpath + '.alg.pkl'):
        print 'loading fasta from pkl file'
        cdna = cPickle.load(open(outpath + '.alg.pkl'))
    else:
        print 'loading from fasta file'
        cdna = SeqGroup(infile, fix_duplicates=False)
        cPickle.dump(cdna, open(outpath + '.alg.pkl', 'wb'), 2)
        print 'loading fasta done'
    print 'start counting'
    spcodons = defaultdict(Counter)
    spbases = defaultdict(Counter)
    spstops = defaultdict(Counter)
    spstarts = defaultdict(Counter)

    for n, (seqid, seq, _) in enumerate(cdna):
        if n % 10000 == 0:
            print "%s sequences processed" % n
        taxid = seqid.split(".")[0]
        # base counts
        spbases[taxid].update(seq)
        # start and stop counts
        start_codon = seq[:3]
        stop_codon = seq[-3:]
        spstarts[taxid].update([start_codon])
        spstops[taxid].update([stop_codon])
        if start_codon in all_starts and stop_codon in all_stops:
            iseq = seq[3:-3]
        elif start_codon in all_starts:
            iseq = seq[3:]
        elif stop_codon in all_stops:
            iseq = seq[:-3]
        else:
            iseq = seq
            
        # iterate over codons
        for i in xrange(0, len(iseq), 3):
            codon = iseq[i: i+3]
            if len(codon) != 3:
                continue
            # codon counts
            spcodons[taxid].update([codon])

    cPickle.dump(spcodons, open(outpath + '.spcodons.pkl', 'wb'), 2)
    cPickle.dump(spbases, open(outpath + '.spbases.pkl', 'wb'), 2)
    cPickle.dump(spstarts, open(outpath + '.spstarts.pkl', 'wb'), 2)
    cPickle.dump(spstops, open(outpath + '.spstops.pkl', 'wb'), 2)

for sp in spbases:
    base_contribution = ambiguity_base(spbases[sp])
    # if there is any ambiguous base we remove ambiguous codons after adding contributions
    if base_contribution:
        ambiguity_codon(spcodons[sp], base_contribution)
        ambiguity_codon(spstarts[sp], base_contribution)
        ambiguity_codon(spstops[sp], base_contribution)

spbasesD = pd.DataFrame.from_dict(spbases, orient='index').fillna(0).loc[tNCBI.get_leaf_names()]
spcodonsD = pd.DataFrame.from_dict(spcodons, orient='index').fillna(0).loc[tNCBI.get_leaf_names(),all_codons]
spstartsD = pd.DataFrame.from_dict(spstarts, orient='index').fillna(0).loc[tNCBI.get_leaf_names(),all_codons]
spstopsD = pd.DataFrame.from_dict(spstops, orient='index').fillna(0).loc[tNCBI.get_leaf_names(),all_codons]


        
