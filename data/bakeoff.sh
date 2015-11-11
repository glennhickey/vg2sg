#!/usr/bin/env bash

#set -e

ASSEMBLY="GRCh38"
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

# get the fasta file associated with a region
function get_region_fa {

	 REGION=$1

	 echo ${REGION}_${ASSEMBLY}/${REGION}.fa
}

for REGION in BRCA1 BRCA2 MHC SMA LRC_KIR
do
	 FA_FILE=`get_region_fa ${REGION}`
	 VCF_FILE=`get_region_vcf ${REGION} 0`
	 COORDS=`get_region_coords ${REGION}`

	 if [[ ! -e ${VCF_FILE} ]] || [[ ! -e ${FA_FILE} ]]
	 then
		  ./fetch1kgpRegion.py ${REGION} --assembly ${ASSEMBLY}
	 fi

	 if [ ! -e "${REGION}_${ASSEMBLY}/${REGION}.vg" ]
	 then
		  vg construct -r ${FA_FILE} -v ${VCF_FILE} -p > ${REGION}_${ASSEMBLY}/${REGION}.vg
	 fi

	 if [ ! -e "${REGION}_${ASSEMBLY}/${REGION}_merge_${WINDOW}.vg" ]
	 then
		  MERGE_VCF_FILE=`get_region_vcf ${REGION} 1`
		  if [ ! -e "${MERGE_VCF_FILE}" ]
		  then
				# merge up blocks of phased snps into their own haplotypes
				gzip -d -c ${VCF_FILE} | vcfgeno2haplo -r ${FA_FILE} -w ${WINDOW} > ${MERGE_VCF_FILE}
		  fi
		  
		  # create compressed vcf
		  bgzip -f ${MERGE_VCF_FILE}
		  tabix -f -p vcf ${MERGE_VCF_FILE}.gz

		  # make vg with -f (keep alleles)
		  vg construct -f -r ${FA_FILE} -v ${VCF_FILE} -p > ${REGION}_${ASSEMBLY}/${REGION}_merge_${WINDOW}.vg
	 fi
done

for REGION in BRCA1 BRCA2 MHC SMA LRC_KIR
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

