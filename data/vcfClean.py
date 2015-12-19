#!/usr/bin/env python2.7

"""
 The 1000 Genomes VCF files have some SNPs which are not congruant with
 the fasta sequence.  This script simply deletes such cases so they 
 don't cause trouble down the road
"""


import argparse, sys, os, os.path, random, subprocess, shutil, itertools
from Bio import SeqIO
from Bio.Seq import Seq

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

def main(args):
    options = parse_args(args)

    in_vcf = open(options.in_vcf, "r")
    in_fa = open(options.in_fa)
    records = list(SeqIO.parse(in_fa, "fasta"))
    assert len(records) == 1
    fa_seq = records[0]
    skip_count = 0
    rec_count = 0

    while True:
        line = in_vcf.readline()
        if line != "":
            skip = False
            if line[0] != "#":
                rec_count += 1
                toks = line.split()
                assert fa_seq.name == toks[0]
                assert len(toks) > 3
                vcf_ref = toks[3]
                start = int(toks[1]) - 1
                fa_ref = ""
                for i in xrange(len(vcf_ref)):
                    fa_ref += str(fa_seq[start + i])
                if fa_ref.upper() != vcf_ref.upper():
                    skip = True
                    skip_count += 1
                    sys.stderr.write("Skipping VCF variant at {} with ref {} because it does not match fasta {}\n".format(toks[1], vcf_ref, fa_ref))
            if not skip:
                sys.stdout.write(line)
        else:
            break

    sys.stderr.write("Skipped {} out of {} records".format(skip_count, rec_count))
    in_vcf.close()
    in_fa.close()
	 
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
