#!/usr/bin/env bash

#
# FROM https://raw.githubusercontent.com/adamnovak/vg/examples/examples/1kgp/1kgp.sh
#


set -e

# This script downloads the 1000 genomes VCFs and the reference genome, and
# turns them into a vg graph. It takes a lot of CPU, memory, network, and disk
# space. It measures just how much CPU and memory it uses with GNU time.

# We're going to download 1000 genomes' VCFs and its reference (GRCh37) in single-chromosome FASTAs.
# They have chromosome FASTAs identified by IDs "1", "2", "X", etc.

# Function to download 1000 genomes VCFs per chromosome
function get_1000g_vcf {
    # This is the contig name we want: "1", "2", "Y"
    CONTIG=$1
    
    # What's the 1000g base url?
    BASE_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
    
    #if [ "$CONTIG" == "Y" ]
    #then
    #    # chrY has a special filename
    #    VCF_URL="${BASE_URL}/ALL.chr${CONTIG}.phase3_integrated_v1a.20130502.genotypes.vcf.gz"
    #else
    #    VCF_URL="${BASE_URL}/ALL.chr${CONTIG}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    #fi
    
    # Actually we will always use this one combined no-samples VCF.
    # Trying to use the with-samples chr1 could use hundreds of gigabytes, because vg loads the whole vcf first.
    VCF_URL="${BASE_URL}/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz"
    
    echo "Retrieving ${VCF_URL}"
    
    # Where do we save it?
    OUTPUT_FILE="vcf/${CONTIG}.vcf.gz"
    
    # We download once and make links
    DOWNLOAD_FILE="vcf/all.vcf.gz"    
    
    # Download the VCF. If any is already downloaded, resume.
    curl -C - -o "${DOWNLOAD_FILE}" "${VCF_URL}"
    # Get the index too
    curl -C - -o "${DOWNLOAD_FILE}.tbi" "${VCF_URL}.tbi"
    
    # Make hard links
    ln "${DOWNLOAD_FILE}" "${OUTPUT_FILE}"
    ln "${DOWNLOAD_FILE}.tbi" "${OUTPUT_FILE}.tbi"
}

# Function to download GRCh37 reference FASTAs per chromosome
function get_GRCh37_fasta {
    # This is the contig name we want: "1", "2", "Y"
    CONTIG=$1
    
    # This is the URL to get it from
    FASTA_URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr${CONTIG}.fa.gz"
    
    # And the place to put it
    OUTPUT_FILE="fa/${CONTIG}.fa"
    OUTPUT_FILE_ZIPPED="${OUTPUT_FILE}.gz"
    
    echo "Retrieving ${FASTA_URL}"
    
    # Go download the zipped version
    # Resume doesn't work properly for this server in http 1.1 if the file is done already.
    curl -C - -0 -o "${OUTPUT_FILE_ZIPPED}" "${FASTA_URL}"

    # We need to strip the "chr" from the record names, and unzip for vg.
    zcat "${OUTPUT_FILE_ZIPPED}" | sed "s/chr${CONTIG}/${CONTIG}/" > ${OUTPUT_FILE}
}

if [ ! -e $0 ]
then
    echo "Please run from the script's directory."
    exit 1
fi

# We want time output on our stdout, but vg progress abr stuff on our stderr
TIME_FILE="$(mktemp)"

# Make the input file directories
mkdir -p vcf
mkdir -p fa

# Make a directory for all the VG files
mkdir -p vg_parts

# Loop over every chromosome
for CONTIG in {1..22} X Y
do

    if [ ! -e "vg_parts/${CONTIG}.vg" ]
    then
    
        # Get the VCF and FASTA
        get_1000g_vcf ${CONTIG}
        get_GRCh37_fasta ${CONTIG}
    
        # Build the vg graph for this contig
        echo "Building VG graph for chromosome ${CONTIG}"
        # We redirect time's time info to a file, and then print it to stdout.
        /usr/bin/time -v -o "${TIME_FILE}" ../../vg construct -r "fa/${CONTIG}.fa" -v "vcf/${CONTIG}.vcf.gz" -R "${CONTIG}" -p > "vg_parts/${CONTIG}.vg"
        cat "${TIME_FILE}"
        
    fi
done

if [ ! -e 1kgp.vg ]
then

    # Now we put them in the same ID space
    echo "Re-numbering nodes to joint ID space..."
    /usr/bin/time -v -o "${TIME_FILE}" ../../vg ids -j vg_parts/*.vg
    cat "${TIME_FILE}"

    # Now we concatenate all the parts into one file
    echo "Collecting subgraphs..."
    cat vg_parts/*.vg > 1kgp.vg
    
fi

# And get some stats
echo "Calculating statistics..."
/usr/bin/time -v -o "${TIME_FILE}" ../../vg stats -z 1kgp.vg
cat "${TIME_FILE}"

rm "${TIME_FILE}"




