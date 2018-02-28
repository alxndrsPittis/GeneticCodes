#!/usr/bin/env python
"""Based on Biopython Cookbook example showing how to parse features from a GenBank file.

This script prints out the mRNA sequences, including stop codons
from GenBank file. It was tested and run with mitochondrion.1.genomic.gbff.gz
but should probably work with any GenBank file.
"""
import sys
from sys import exit
import gzip
from Bio import GenBank
#from Bio.Seq import MutableSeq
#from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from ete3 import NCBITaxa

CorrectDict = {"Leptocephalus sp. 'type II larva' (Smith, 1989)": "Leptocephalus sp. 'type II larva'",
               "Humphaplotropis culaishanensis" : "Humphaplotropis culaishanensis (nomen nudum)",
               "Paraglypturus tonganus" : "Paraglypturus tonganus (nomen nudum)",
               "Hoploplana elisabelloi" : "Hoploplana elisabelloi (nomen nudum)",
               "Palpitomonas bilix Eukaryota." : "Palpitomonas bilix",
               "Eukaryota sp. BB2 Eukaryota." : "Eukaryota sp. BB2",
               "Ancoracysta twista Eukaryota." : "Ancoracysta twista"
    }

# --- load a parser and iterator for our GenBank file
gb_handle = gzip.open(sys.argv[1], "r")
# -- a parser that will give you back SeqFeature objects
feature_parser = GenBank.FeatureParser()
iterator = GenBank.Iterator(gb_handle, feature_parser)

# load taxonomy for taxids
ncbi = NCBITaxa()

# output using prefix
out_1 = open("%s.cdna.fasta" % sys.argv[2], "w")
out_2 = open("%s.codons.tab" % sys.argv[2], "w")

strands = []
stop_codons = []
prot_ids = []
excluded = []
missing_id = []

# begin iterating through the file and getting GenBank records
while 1:
    # get a SeqFeature object for the next GenBank record. When we run
    # out of records in the file, cur_entry will be None
    cur_entry = iterator.next()

    if cur_entry is None:
        break

    nid = cur_entry.id
    
    organism = cur_entry.annotations['organism']
    # name exception..
    try:
        taxid = ncbi.get_name_translator([organism])[organism][0]
    except KeyError:
        if organism in CorrectDict:
            taxid = ncbi.get_name_translator([CorrectDict[organism]])[CorrectDict[organism]][0]
        else:
            print cur_entry.annotations
            try:
                correct_organism = raw_input("name for %s\n" % organism)
                taxid = ncbi.get_name_translator([correct_organism])[correct_organism][0]
            except KeyError:
                taxid = int(raw_input("taxid"))
    #print "Printing cDNA info for %s" % nid
    # loop through all of the features for the entry
    for feature in cur_entry.features:
        # when we've got CDS features, parse the info out of them
        if feature.type == "CDS":
            try:
                prot_id = feature.qualifiers['protein_id'][0]
            except KeyError:
                locus_tag = feature.qualifiers['locus_tag'][0]
                print "%s missing id" % locus_tag
                missing_id.append(locus_tag)
                continue
            seqid = "%s.%s" % (taxid, prot_id)
            trail_codon = False
            if 'transl_table' in feature.qualifiers:
                transl_table_id = int(feature.qualifiers['transl_table'][0])
            else:
                transl_table_id = 1

            # Errors, wrong translations etc.
            if prot_id in ['YP_009370766.1']:
                continue
                
            #transl_table = CodonTable.unambiguous_dna_by_id[transl_table_id]
            # stop_codons = transl_table.stop_codons
            # stop_codons = ['TAA', 'TAG', 'AGA']
            try:
                cdna = feature.location.extract(cur_entry.seq)
            except ValueError:
                print "%s excluded" % prot_id
                excluded.append(prot_id)
                continue
            # some cds' translation starts in 2 or 3 position
            if 'codon_start' in feature.qualifiers:
                codon_start = int(feature.qualifiers['codon_start'][0])
                if codon_start != 1:
                    cdna = cdna[(codon_start-1):]

            # keep them to modify if error in translation above
            location_start = feature.location.start
            location_end = feature.location.end
                    
            # case where one or more nucl missing while last codon(s) translated
            if len(cdna) / 3. < len(feature.qualifiers['translation'][0]):
                diff = (len(feature.qualifiers['translation'][0]) * 3) - len(cdna)
                if feature.location.strand == 1:
                    # add what's missing from end
                    location_end += diff
                    if location_end >= len(cur_entry.seq):
                        missing_nucl = cur_entry.seq[feature.location.end:] + \
                                       cur_entry.seq[:(location_end-len(cur_entry.seq))]
                        # New end location, for correct stop codon indices
                        location_end = location_end-len(cur_entry.seq)
                    else:
                        missing_nucl = cur_entry.seq[feature.location.end:location_end]
                    cdna = cdna + missing_nucl
                else:
                    # subtract from start
                    location_start -= diff
                    if location_start < 0:
                        # take it from opposite end
                        missing_nucl = cur_entry.seq[location_start:feature.location.start].reverse_complement() + \
                                       cur_entry.seq[location_start:].reverse_complement()
                        location_start = len(cur_entry.seq) + location_start
                    else:
                        missing_nucl = cur_entry.seq[location_start:feature.location.start].reverse_complement()
                    cdna = cdna + missing_nucl
                
            # Get the trailing nucleotides, not in any codon, usually A's added in mRNA
            # Nucleotides in cdna not a multiple of 3
            remainder = len(cdna) % 3

            if len(cdna) / 3. == len(feature.qualifiers['translation'][0]) + 1:
                # Remainder == 0 if condition is met
                # When reverse gene starts at 0, stop codon is included
                full_cdna = cdna
                stop_codon = cdna[-3:]
            elif len(cdna) > (len(feature.qualifiers['translation'][0]) * 3) + 3:
                # if equal, previous condition would be met.
                # cDNA has excessive nuclotides, removed according to translation..
                nexcessive = (len(cdna) - len(feature.qualifiers['translation'][0]) * 3)
                full_cdna = cdna[:(-nexcessive + 3)]
                stop_codon = full_cdna[-3:]
            else:
                # Stop codon not included or remainder present
                # if remainder and "note" not in feature.qualifiers and len(cur_entry.seq) != location_end:
                #     exit("no note")
                strand = feature.location.strand
                if "note" in feature.qualifiers and \
                   "TAA stop codon is completed by the addition of 3' A residues to the mRNA" in str(feature.qualifiers['note']):
                    if not remainder:
                        if strand == 1:
                            stop_codon = cur_entry.seq[location_end:location_end+1] + "AA"
                            trail_codon = cur_entry.seq[location_end:location_end+3]
                            full_cdna = cdna + stop_codon
                            #full_cdna = cur_entry.seq[location_start:location_end+1] + "AA"
                            #stop_codon = full_cdna[-3:]
                            #trail_codon = cur_entry.seq[location_end:location_end+3]
                        else:
                            stop_codon = cur_entry.seq[location_start-1:location_start].reverse_complement() + "AA"
                            trail_codon = cur_entry.seq[location_start-3:location_start].reverse_complement()
                            full_cdna = cdna + stop_codon                            
                            #full_cdna = cur_entry.seq[location_start-1:location_end].reverse_complement() + "AA"
                            #stop_codon = full_cdna[-3:]
                            #trail_codon = cur_entry.seq[location_start-3:location_start].reverse_complement()
                        if not stop_codon.startswith("T"):
                            exit("trail nucl starts not with T")
                    else:
                        stop_codon = cdna[-remainder:] + "A" * (3-remainder)
                        full_cdna = cdna[:-remainder] + stop_codon
                        # keep codons with trailing nucl also..
                        if strand == -1:
                            start_pos = location_start-3+remainder
                            end_pos = location_start+remainder
                            # account for cases where inclomplete stop at 0
                            if start_pos < 0:
                                trail_codon = cur_entry.seq[0:end_pos].reverse_complement() +\
                                              cur_entry.seq[start_pos:].reverse_complement()
                            else:
                                trail_codon = cur_entry.seq[start_pos:end_pos].reverse_complement()
                        else:
                            start_pos = location_end-remainder
                            end_pos = location_end-remainder+3
                            trail_codon = cur_entry.seq[start_pos:end_pos]
                            # account for cases where inclomplete stop at end
                            # we take positions from start, as genome circular
                            if len(trail_codon) != 3 and strand == 1:
                                trail_codon += cur_entry.seq[:(3-len(trail_codon))]
                        if len(trail_codon) != 3:
                            exit("inclomplete trailing stop codon")
                else:
                    if strand == -1:
                        # - strand
                        start_pos = location_start-3+remainder
                        end_pos = location_start+remainder
                        # account for cases where inclomplete stop at 0
                        if start_pos < 0:
                            stop_codon = cur_entry.seq[0:end_pos].reverse_complement() +\
                                          cur_entry.seq[start_pos:].reverse_complement()
                        else:
                            stop_codon = cur_entry.seq[start_pos:end_pos].reverse_complement()
                    else:
                        # + strand
                        start_pos = location_end-remainder
                        end_pos = location_end-remainder+3
                        stop_codon = cur_entry.seq[start_pos:end_pos]
                        # if at the end, take positions from start, as genome circular
                        if len(stop_codon) != 3:
                            stop_codon += cur_entry.seq[:(3-len(stop_codon))]
                    if remainder:
                        full_cdna = cdna[:-remainder] + stop_codon
                    else:
                        full_cdna = cdna + stop_codon
                        
            if len(stop_codon) != 3:
                exit("incomplete or long stop codon")
            if not len(full_cdna) / 3. == len(feature.qualifiers['translation'][0]) + 1:
                print prot_id
                print feature.location.extract(cur_entry.seq).translate()
                print feature.qualifiers['translation'][0]
                exit("cdna length issue")
            # exclude NNN sequences
            if len(set(full_cdna)) == 1:
                continue


            stop_codons.append(stop_codon)
            # cdna fasta file
            print >>out_1, ">%s" % seqid
            print >>out_1, full_cdna
            # codon dictionary
            if trail_codon:
                print >>out_2, "%s\t%s\t%s\t%s" % (seqid, transl_table_id, stop_codon, trail_codon)
            else:
                print >>out_2, "%s\t%s\t%s" % (seqid, transl_table_id, stop_codon)

if excluded:
    print "%s excluded entries:" % len(excluded)
    print ",".join(excluded)

if missing_id:
    print "%s missing id entries:" % len(missing_id)
    print ",".join(missing_id)
    
out_1.close()
out_2.close()
