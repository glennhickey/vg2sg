#!/usr/bin/env python2.7

# based on Adam Novak's excellent scripts for doing similar stuff
# (fetchRegion.py, 1kgp.sh: mirrored in this directory)


"""
fetch1kgpRegion.py: Fetch the sequence FASTA data and 1000 genomes VCF file for
region (like "LRC_KIR" or "MHC") by name.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools
import collections, urllib2, shutil, subprocess, glob, doctest

from Bio import AlignIO, SeqIO, Align, Entrez

import tsv

from fetchRegion import get_region_info, get_ucsc_name


def parse_args(args):
    """
    Takes in the command-line arguments list (args), and returns a nice argparse
    result with fields for all the options.
    
    Borrows heavily from the argparse documentation examples:
    <http://docs.python.org/library/argparse.html>
    """
    
    # Construct the parser (which is stored in parser)
    # Module docstring lives in __doc__
    # See http://python-forum.com/pythonforum/viewtopic.php?f=3&t=36847
    # And a formatter class so our examples in the docstring look good. Isn't it
    # convenient how we already wrapped it to 80 characters?
    # See http://docs.python.org/library/argparse.html#formatter-class
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # General options
    parser.add_argument("region",
        help="name of the region to download, and the output directory")
    parser.add_argument("--assembly_url", 
        default=("ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/"
                 "vertebrates_mammals/Homo_sapiens/GRCh37.p13/"),
                 
        # default is GRC37.  This would be dir for 38:
        #default=("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/"
        #    "GCA_000001405.17_GRCh38.p2/"
        #    "GCA_000001405.17_GRCh38.p2_assembly_structure"),

        help="URL for the assembly, containing genomic_region_definitions.txt")
    parser.add_argument("--email", default="soe@soe.ucsc.edu",
        help="E-mail address to report to Entrez")

    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

# Function to download 1000 genomes VCFs per chromosome
def get_1000g_vcf(CONTIG, SLICE_PATH, BED_PATH):

    if not os.path.exists("vcf"):
        os.makedirs("vcf")
        
    # What's the 1000g base url?
    BASE_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
        
    # Actually we will always use this one combined no-samples VCF.
    # Trying to use the with-samples chr1 could use hundreds of gigabytes, because vg loads the whole vcf first.
    VCF_URL="{}/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz".format(BASE_URL)
    
    print "Retrieving {}".format(VCF_URL)
    
    # Where do we save it?
    OUTPUT_FILE="vcf/{}.vcf.gz".format(CONTIG)
    
    # We download once and make links
    DOWNLOAD_FILE="vcf/all.vcf.gz"
    VCF_FILE="vcf/all.vcf"
    
    # Download the VCF. If any is already downloaded, resume.
    os.system("curl -C - -o {} {}".format(DOWNLOAD_FILE, VCF_URL))
    # Get the index too
    os.system("curl -C - -o {}.tbi {}.tbi".format(DOWNLOAD_FILE, VCF_URL))
    
    # Make hard links
    os.system("ln {} {}".format(DOWNLOAD_FILE, OUTPUT_FILE))
    os.system("ln {}.tbi {}.tbi".format(DOWNLOAD_FILE, OUTPUT_FILE))

    # Uncompress vcf for Bedtools
    os.system("gzip -d {}".format(DOWNLOAD_FILE))

    # Intersect region
    os.system("intersectBed -a {} -b {} > {}".format(BED_PATH, VCF_FILE,
                                                     SLICE_PATH))

    # Keep vcf compressed
    os.system("gzip {}".format(VCF_FILE))


# Function to download GRCh37 reference FASTAs per chromosome
def get_GRCh37_fasta(CONTIG, SLICE_PATH, BED_PATH):

    if not os.path.exists("fa"):
        os.makedirs("fa")

    # This is the URL to get it from
    FASTA_URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr{}.fa.gz".format(CONTIG)
    
    # And the place to put it
    OUTPUT_FILE="fa/{}.fa".format(CONTIG)
    OUTPUT_FILE_ZIPPED="{}.gz".format(OUTPUT_FILE)
    
    print "Retrieving {}".format(FASTA_URL)
    
    # Go download the zipped version
    # Resume doesn't work properly for this server in http 1.1 if the file is done already.
    os.system("curl -C - -0 -o {} {}".format(OUTPUT_FILE_ZIPPED, FASTA_URL))
              
    # We need to strip the "chr" from the record names, and unzip for vg.
    os.system("cat {} | zcat | sed \"s/chr{}/{}/\" > {}".format(OUTPUT_FILE_ZIPPED, CONTIG, CONTIG, OUTPUT_FILE))

    # now get a FASTA file for just the contig
    print "bedtools getfasta -fi {} -bed {} -fo {}".format(OUTPUT_FILE, BED_PATH, SLICE_PATH)
    os.system("bedtools getfasta -fi {} -bed {} -fo {}".format(OUTPUT_FILE, BED_PATH, SLICE_PATH))
    
def main(args):
    
    options = parse_args(args) # This holds the nicely-parsed options object

    # Set Entrez e-mail
    Entrez.email = options.email

    # Go get the region of the reference we're talking about. Starts and ends
    # are 1-based.
    ref_acc, ref_start, ref_end = get_region_info(options.region,
        options.assembly_url)

    # get the ucsc chromosome name
    ref_chr = get_ucsc_name(ref_acc)

    # use contig names that don't have chr throughout
    contig = ref_chr.replace("chr", "")
        
    # Make our output directory
    if not os.path.exists(options.region):
        os.makedirs(options.region)

    # Make our BED region
    bed_path = os.path.join(options.region, options.region + ".bed")
    bed_file = open(bed_path, "w")
    bed_file.write("{}\t{}\t{}\t\n".format(contig, ref_start, ref_end))

    # Download 1000 Genomes data and slice region out
    fa_path = os.path.join(options.region, options.region + ".fa")
    vcf_path = os.path.join(options.region, options.region + ".vcf")
    get_GRCh37_fasta(contig, fa_path, bed_path)
    get_1000g_vcf(contig, vcf_path, bed_path)
                    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
