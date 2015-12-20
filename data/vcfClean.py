#!/usr/bin/env python2.7

"""
 The 1000 Genomes VCF files have some SNPs which are reversed 
 in LRC_KIR.  ie the reference in the VCF is the reverse complement
 of teh FASTA.  checking a few in dbsnp it seems that the ALTs are
 reversed as well, so the snp can be corrected in the VCF by reversing
 everything. 
"""


import argparse, sys, os, os.path, random, subprocess, shutil, itertools
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("in_vcf", type=str,
                        help="Input vcf file"),
    parser.add_argument("in_fa", type=str,
                        help="Input fa")
    
    args = args[1:]
    options = parser.parse_args(args)
    return options

def flip_record(toks):
    vcf_ref = toks[3]
    vcf_alts = toks[4].split(",")
    fix_ref = reverse_complement(vcf_ref)
    fix_alts = ",".join([reverse_complement(x) for x in vcf_alts])
    fix_line = toks[0:3] + [fix_ref, fix_alts] + toks[5:]
    fix_line = "\t".join(fix_line)
    return fix_line + "\n" 
    
def main(args):
    options = parse_args(args)

    in_vcf = open(options.in_vcf, "r")
    in_fa = open(options.in_fa)
    records = list(SeqIO.parse(in_fa, "fasta"))
    assert len(records) == 1
    fa_seq = records[0]
    skip_count = 0
    rec_count = 0
    flip_count = 0

    while True:
        line = in_vcf.readline()
        if line == "":
            break
        elif line[0] == "#":
            sys.stdout.write(line)
        else:
            rec_count += 1
            toks = line.split()
            assert fa_seq.name == toks[0]
            assert len(toks) > 3
            vcf_ref = toks[3]
            vcf_alts = toks[4]
            start = int(toks[1]) - 1
            fa_ref = ""
            for i in xrange(len(vcf_ref)):
                fa_ref += str(fa_seq[start + i])
            if fa_ref.upper() != vcf_ref.upper():
                flip_count += 1
                line = flip_record(toks)
                vcf_ref = line.split()[3]
            if fa_ref.upper() != vcf_ref.upper():
                skip_count += 1
                sys.stderr.write("Skipping VCF variant at {} with ref {} because it does not match fasta {} and is not fixable by reversing\n".format(toks[1], toks[3], fa_ref))
            else:
                sys.stdout.write(line)

    sys.stderr.write("Skipped {} and flipped {} out of {} records".format(skip_count, flip_count, rec_count))
    in_vcf.close()
    in_fa.close()
	 
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
