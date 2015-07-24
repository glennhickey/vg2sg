#!/usr/bin/env bash

#set -e

VCF="vcf/all.vcf.gz"

# get the fasta file associated with a region
function get_region_fa {

	 REGION=$1

	 # get the bed path
	 BED=${REGION}/${REGION}.bed

	 # get contig
	 CONTIG=`cat ${BED} | awk '{print $1}'`

	 echo "fa/${CONTIG}.fa"
}

# get the vcf -R coordinates for a region
function get_region_coords {

	 REGION=$1

	 # get the bed path
	 BED=${REGION}/${REGION}.bed

	 # get contig
	 CONTIG=`cat ${BED} | awk '{print $1}'`

	 #get coordinates
	 START=`cat ${BED} | awk '{print $2}'`

	 #note VCF seems to want inclusive end
	 END=`cat ${BED} | awk '{print $3-1}'`

	 echo "${CONTIG}:${START}-${END}"
}

for REGION in MHC SMA
do
	 if [ ! -e "${REGION}/${REGION}.vg" ]
	 then
		  ./fetch1kgpRegion.py ${REGION}
		  FA_FILE=`get_region_fa ${REGION}`
		  COORDS=`get_region_coords ${REGION}`
		  vg construct -r ${FA_FILE} -R ${COORDS} -v ${VCF} -p > ${REGION}/${REGION}.vg
	 fi
done

for REGION in MHC SMA
do
	 if [ ! -e "${REGION}/database.sql" ]
	 then
		  ../vg2sg ${REGION}/${REGION}.vg ${REGION}/database.fa ${REGION}/database.sql -s
	 fi
done

