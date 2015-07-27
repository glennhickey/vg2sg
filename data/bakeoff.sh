#!/usr/bin/env bash

#set -e

ASSEMBLY="GRCh38"

#get the vcf file associated with a region
function get_region_vcf {

	 REGION=$1

	 # get the bed path
	 BED=${REGION}_${ASSEMBLY}/${REGION}.bed

	 # get contig
	 CONTIG=`cat ${BED} | awk '{print $1}'`

	 echo "vcf_${ASSEMBLY}/${CONTIG}.vcf.gz"
}

# get the fasta file associated with a region
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

for REGION in MHC SMA LRC_KIR
do
	 if [ ! -e "${REGION}_${ASSEMBLY}/${REGION}.vg" ]
	 then
		  ./fetch1kgpRegion.py ${REGION} --assembly ${ASSEMBLY}
		  FA_FILE=`get_region_fa ${REGION}`
		  VCF_FILE=`get_region_vcf ${REGION}`
		  COORDS=`get_region_coords ${REGION}`
		  vg construct -r ${FA_FILE} -R ${COORDS} -v ${VCF_FILE} -p > ${REGION}_${ASSEMBLY}/${REGION}.vg
	 fi
done

for REGION in MHC SMA LRC_KIR
do
	 if [ ! -e "${REGION}_${ASSEMBLY}/database.sql" ]
	 then
		  vg2sg ${REGION}_${ASSEMBLY}/${REGION}.vg ${REGION}_${ASSEMBLY}/database.fa ${REGION}_${ASSEMBLY}/database.sql -s
	 fi
done

