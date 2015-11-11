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
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("region",
            help="name of the region to download, and the output directory")
    parser.add_argument("--assembly", default="GRCh38",
                        help="Either GRCh37 or GRCh38")
    parser.add_argument("--email", default="soe@soe.ucsc.edu",
                        help="E-mail address to report to Entrez")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]

    options = parser.parse_args(args)
    
    # options.assembly_url = URL for the assembly, containing
    # genomic_region_definitions.txt")
    if options.assembly == "GRCh38":
        options.assembly_url = ("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/"
                "GCA_000001405.17_GRCh38.p2/"
                "GCA_000001405.17_GRCh38.p2_assembly_structure")
    elif options.assembly == "GRCh37":
        options.assembly_url = ("ftp://ftp.ncbi.nlm.nih.gov/genbank/"
                                "genomes/Eukaryotes/"
                                "vertebrates_mammals/Homo_sapiens/GRCh37.p13/")
    else:
        raise "Invalid argument for --assembly option"
    
    return options

# Function to download 1000 genomes VCFs per chromosome
def get_1000g_vcf(CONTIG, ASSEMBLY):

    if not os.path.exists("vcf_{}".format(ASSEMBLY)):
        os.makedirs("vcf_{}".format(ASSEMBLY))

    # What's the 1000g base url?
    BASE_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
        
    if ASSEMBLY == "GRCh37":
        if CONTIG == "Y":
            # chrY has a special filename
            VCF_URL=("{}/ALL.chr{}.phase3_integrated_v1a.20130502."
                "genotypes.vcf.gz").format(BASE_URL, CONTIG)
        else:
            VCF_URL=("{}/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a."
                     "20130502.genotypes.vcf.gz").format(BASE_URL, CONTIG)

    else:
        assert ASSEMBLY == "GRCh38"
        assert CONTIG != "X" and CONTIG != "Y"
        BASE_URL += "/supporting/GRCh38_positions"
        VCF_URL = ("{}/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v3plus"
                   "_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz"
                   ).format(BASE_URL, CONTIG)
        
    print "Retrieving {}".format(VCF_URL)
        
    # We download once and make links
    DOWNLOAD_FILE="vcf_{}/{}.vcf.gz".format(ASSEMBLY, CONTIG)
    
    # Download the VCF. If any is already downloaded, resume.
    os.system("curl -C - -o {} {}".format(DOWNLOAD_FILE, VCF_URL))
    # Get the index too
    os.system("curl -C - -o {}.tbi {}.tbi".format(DOWNLOAD_FILE, VCF_URL))

 # Function to download reference FASTAs per chromosome
def get_fasta(CONTIG, ASSEMBLY):

    if not os.path.exists("fa_{}".format(ASSEMBLY)):
        os.makedirs("fa_{}".format(ASSEMBLY))

    if ASSEMBLY == "GRCh38":
        NAME = "hg38"
    else:
        assert ASSEMBLY == "GRCh37"
        NAME = "hg19"
        
    # This is the URL to get it from
    FASTA_URL=("http://hgdownload.cse.ucsc.edu/goldenPath/{}/"
               "chromosomes/chr{}.fa.gz").format(NAME, CONTIG)
    
    # And the place to put it
    OUTPUT_FILE="fa_{}/{}.fa".format(ASSEMBLY, CONTIG)
    OUTPUT_FILE_ZIPPED="{}.gz".format(OUTPUT_FILE)
    
    print "Retrieving {}".format(FASTA_URL)
    
    # Go download the zipped version
    # Resume doesn't work properly for this server in http 1.1 if the file is done already.
    os.system("curl -C - -0 -o {} {}".format(OUTPUT_FILE_ZIPPED, FASTA_URL))
              
    # We need to strip the "chr" from the record names, and unzip for vg.
    os.system("cat {} | zcat | sed \"s/chr{}/{}/\" > {}".format(
        OUTPUT_FILE_ZIPPED, CONTIG, CONTIG, OUTPUT_FILE))

# return hardcoded brca coordinates (get_region_info() doesn't know brca)
def get_brca_info(region, assembly):
    # brca1 http://www.ncbi.nlm.nih.gov/gene/672
    # brca2 http://www.ncbi.nlm.nih.gov/gene/675
    if assembly == "GRCh38":
        if region == "BRCA1":
            return "CM000679.2", 43044295, 43125483
        elif region == "BRCA2":
            return "CM000675.2", 32314862, 32399850

    elif assembly == "GRCh37":
        if region == "BRCA1":
            return "CM000679.1", 41196312, 41277500
        elif region == "BRCA2":
            return "CM000675.1", 32889617, 32973809

    assert False

def main(args):
    
    options = parse_args(args) # This holds the nicely-parsed options object

    # Set Entrez e-mail
    Entrez.email = options.email

    # Go get the region of the reference we're talking about. Starts and ends
    # are 1-based.
    if options.region == "BRCA1" or options.region == "BRCA2":
        ref_acc, ref_start, ref_end = get_brca_info(options.region,                                                    options.assembly)
    else:
        ref_acc, ref_start, ref_end = get_region_info(options.region,
            options.assembly_url)

    # get the ucsc chromosome name
    ref_chr = get_ucsc_name(ref_acc)

    # use contig names that don't have chr throughout
    contig = ref_chr.replace("chr", "")
        
    # Make our output directory
    reg_dir = "{}_{}".format(options.region, options.assembly)
    if not os.path.exists(reg_dir):
        os.makedirs(reg_dir)

    # Make our BED region
    bed_path = os.path.join(reg_dir, options.region + ".bed")
    bed_file = open(bed_path, "w")
    # Note: we keep BED coordinates 0-based:
    bed_file.write("{}\t{}\t{}\t\n".format(contig, ref_start - 1, ref_end))

    # Download 1000 Genomes data and slice region out
    get_fasta(contig, options.assembly)
    get_1000g_vcf(contig, options.assembly)
                    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
