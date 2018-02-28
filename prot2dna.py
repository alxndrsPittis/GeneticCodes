#!/usr/bin/env python
"""
Codon aware reverse translation of protein alignment to dna.
Assumes * is not present in protein alns

"""
import sys
import os
from sys import exit
from glob import glob
from ete3 import SeqGroup, parser

path = sys.argv[1] + "*clustalo"
infiles = glob(path)
print "%s infiles" % len(infiles)

F = parser.fasta.read_fasta(sys.argv[2])

for infile in infiles:
    print infile
    if os.stat(infile).st_size == 0:
        continue
    alg_aa = SeqGroup(infile)
    alg_dna = SeqGroup()

    for name, seq, _ in alg_aa:
        try:
            cdna = F.id2seq[F.name2id[name]]
        except KeyError:
            print "cdna for %s not found" % name
            continue
        cdna_aln = ""
        for pos in seq:
            if pos != "-":
                cdna_aln += cdna[:3]
                cdna = cdna[3:]
            else:
                cdna_aln += "---"
        # Last the stop codon
        cdna_aln += cdna[:3]
        alg_dna.set_seq(name, cdna_aln)
    print "Input protein alignment contains %s aa sequences" % len(alg_aa)
    print "Output cdna alignment contains %s cdna sequences" % len(alg_dna)
    print 
    alg_dna.write(outfile=infile.replace(".clustalo", ".clustalo.cdna.aln"))









                
