#!/usr/bin/env bash

#set -e

ASSEMBLY="GRCh38"
REGIONS=( "BRCA1" "BRCA2" "SMA" "LRC_KIR" "MHC" )
WINDOW=100

#get the vcf file associated with a region
function get_region_vcf {

	 REGION=$1
	 MERGE=$2

	 if [ $MERGE != "1" ]
	 then
		  echo "${REGION}_${ASSEMBLY}/${REGION}.vcf.gz"
	 else
		  echo "${REGION}_${ASSEMBLY}/${REGION}_merge_${WINDOW}.vcf.gz"
	 fi
}

# get the fasta file associated with a region's contig
function get_region_fa {

    REGION=$1

    # get the bed path
    BED=${REGION}_${ASSEMBLY}/${REGION}.bed

    # get contig
    CONTIG=`cat ${BED} | awk '{print $1}'`

    echo "fa_${ASSEMBLY}/${CONTIG}.fa"
}

# get the vcf -R coordinates for a region
function get_region_coords {

    REGION=$1

    # get the bed path
    BED=${REGION}_${ASSEMBLY}/${REGION}.bed

    # get contig
    CONTIG=`cat ${BED} | awk '{print $1}'`

    #get coordinates (convert BED into 1-based, inclusive)
    START=`cat ${BED} | awk '{print $2+1}'`

    #note VCF seems to want inclusive end
    END=`cat ${BED} | awk '{print $3}'`

    echo "${CONTIG}:${START}-${END}"
}

for REGION in "${REGIONS[@]}"
do
	 ./fetch1kgpRegion.py ${REGION} --assembly ${ASSEMBLY}
	 FA_FILE=`get_region_fa ${REGION}`
	 VCF_FILE=`get_region_vcf ${REGION} 0`
	 MERGE_VCF_FILE=`get_region_vcf ${REGION} 1`
	 COORDS=`get_region_coords ${REGION}`

	 if [ ! -e "${REGION}_${ASSEMBLY}/${REGION}.vg" ]
	 then
		  vg construct -r ${FA_FILE} -v ${VCF_FILE} -p > ${REGION}_${ASSEMBLY}/${REGION}.vg
	 fi

	 if [ ! -e "${REGION}_${ASSEMBLY}/${REGION}_merge_${WINDOW}.vg" ]
	 then
		  # merge up blocks of phased snps into their own haplotypes
		  gzip -d -c ${VCF_FILE} | vcfgeno2haplo -r ${FA_FILE} -w ${WINDOW} > ${MERGE_VCF_FILE%??}
		  
		  # create compressed vcf
		  bgzip -f ${MERGE_VCF_FILE%??}
		  tabix -f -p vcf ${MERGE_VCF_FILE}
		  
		  # make vg with -f (keep alleles)
		  vg construct -f -r ${FA_FILE} -v ${MERGE_VCF_FILE} -p > ${REGION}_${ASSEMBLY}/${REGION}_merge_${WINDOW}.vg
	 fi
done

for REGION in "${REGIONS[@]}"
do
	 if [ ! -e "${REGION}_${ASSEMBLY}/database.sql" ]
	 then
		  vg2sg ${REGION}_${ASSEMBLY}/${REGION}.vg ${REGION}_${ASSEMBLY}/database.fa ${REGION}_${ASSEMBLY}/database.sql -s
	 fi
	 if [ ! -e "${REGION}_${ASSEMBLY}/database_merge_${WINDOW}.sql" ]
	 then
		  vg2sg ${REGION}_${ASSEMBLY}/${REGION}_merge_${WINDOW}.vg ${REGION}_${ASSEMBLY}/database_merge_${WINDOW}.fa ${REGION}_${ASSEMBLY}/database_merge_${WINDOW}.sql -s
	 fi

done

