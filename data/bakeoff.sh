#!/usr/bin/env bash

#set -e

ASSEMBLY="GRCh38"
REGIONS=( "BRCA1" "BRCA2" "SMA" "LRC_KIR" "MHC" )
#REGIONS=( "BRCA1" "BRCA2" )
WINDOWS=( 30 50 )
KMER=27
EDGE=3
# from ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/analysis/Illumina_PlatinumGenomes_NA12877_NA12878_09162015/IlluminaPlatinumGenomes-user-guide.pdf
PEDI_SAMPLES="NA12889 NA12890 NA12891 NA12892 NA12877 NA12878 NA12879 NA12880 NA12881 NA12882 NA12883 NA12884 NA12885 NA12886 NA12887 NA12888 NA12893"

#get the vcf file associated with a region
function get_region_vcf {

	 REGION=$1

	 # get the bed path
    BED=${REGION}_${ASSEMBLY}/${REGION}.bed

    # get contig
    CONTIG=`cat ${BED} | awk '{print $1}'`

	 echo "vcf_${ASSEMBLY}/${CONTIG}.vcf.gz"
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

# like above but just the start
function get_start_coord {

	 REGION=$1

	 # get the bed path
	 BED=${REGION}_${ASSEMBLY}/${REGION}.bed

	 #get coordinates (convert BED into 1-based, inclusive)
	 START=`cat ${BED} | awk '{print $2+1}'`

	 echo "${START}"
}

# like above but just the chrom
function get_chrom_coord {

	 REGION=$1

	 # get the bed path
	 BED=${REGION}_${ASSEMBLY}/${REGION}.bed

	 #get coordinates (convert BED into 1-based, inclusive)
	 CHROM=`cat ${BED} | awk '{print $1}'`

	 echo "${CHROM}"
}


for REGION in "${REGIONS[@]}"
do
	 ./fetch1kgpRegion.py ${REGION} --assembly ${ASSEMBLY}
	 FA_FILE=`get_region_fa ${REGION}`
	 VCF_FILE=`get_region_vcf ${REGION} 0`
	 COORDS=`get_region_coords ${REGION}`
	 START_COORD=`get_start_coord ${REGION}`
	 REGION_VCF=${REGION}_${ASSEMBLY}/${REGION}.vcf

	 # slice the region out of the vcf
	 bcftools view ${VCF_FILE} -r ${COORDS} > ${REGION_VCF}.raw
	 # remove pedigree samples
	 vcfremovesamples ${REGION_VCF}.raw ${PEDI_SAMPLES} | vcffixup - | vcffilter -f 'AC > 0' > ${REGION_VCF}.ped
	 # fix some bugs in 1000 genomes liftover where strands get bungled
	 ./vcfClean.py ${REGION_VCF}.ped ${FA_FILE} > ${REGION_VCF} 2> ${REGION}_${ASSEMBLY}/${REGION}.vcfclean.log
	 bgzip ${REGION_VCF} -c > ${REGION_VCF}.gz
	 tabix -f -p vcf ${REGION_VCF}.gz

	 # make a normal vg
	 echo "vg construct -r ${FA_FILE} -v ${REGION_VCF}.gz -R ${COORDS} -p > ${REGION}_${ASSEMBLY}/${REGION}.vg"
	 vg construct -r ${FA_FILE} -v ${REGION_VCF}.gz -R ${COORDS} -p > ${REGION}_${ASSEMBLY}/${REGION}.vg
	 vg index -s -k ${KMER} -e ${EDGE} ${REGION}_${ASSEMBLY}/${REGION}.vg

	 # make a flat vg for snpBridge input
	 vg construct -f -r ${FA_FILE} -v ${REGION_VCF}.gz -R ${COORDS} -p > ${REGION}_${ASSEMBLY}/${REGION}_f.vg
	 vg index -s -k ${KMER} -e ${EDGE} ${REGION}_${ASSEMBLY}/${REGION}_f.vg

	 # make a bridged vg for each window
	 for WINDOW in "${WINDOWS[@]}"
	 do
		  snpBridge ${REGION}_${ASSEMBLY}/${REGION}_f.vg ${REGION_VCF} -o ${START_COORD} -w ${WINDOW} > ${REGION}_${ASSEMBLY}/${REGION}_bridge_${WINDOW}.vg 2> ${REGION}_${ASSEMBLY}/${REGION}_bridge_${WINDOW}.log
		  vg index -s -k ${KMER} -e ${EDGE} ${REGION}_${ASSEMBLY}/${REGION}_bridge_${WINDOW}.vg
		  vg compare ${REGION}_${ASSEMBLY}/${REGION}.vg  ${REGION}_${ASSEMBLY}/${REGION}_bridge_${WINDOW}.vg > ${REGION}_${ASSEMBLY}/${REGION}_bridge_${WINDOW}_comp.txt
	 done

done

for REGION in "${REGIONS[@]}"
do
	 CHROM=`get_chrom_coord ${REGION}`

	 # rename the path to "ref"
	 mv ${REGION}_${ASSEMBLY}/${REGION}.vg ${REGION}_${ASSEMBLY}/${REGION}_chrom_name.vg
	 ./vgRenamePath.sh ${REGION}_${ASSEMBLY}/${REGION}_chrom_name.vg ${CHROM} ref > ${REGION}_${ASSEMBLY}/${REGION}.vg

	 # to side graph sql (-s option because paths dont cover snps)
	 vg2sg ${REGION}_${ASSEMBLY}/${REGION}.vg ${REGION}_${ASSEMBLY}/database.fa ${REGION}_${ASSEMBLY}/database.sql -s

	 for WINDOW in "${WINDOWS[@]}"
	 do
		  # rename the path to "ref"
		  mv ${REGION}_${ASSEMBLY}/${REGION}_bridge_${WINDOW}.vg ${REGION}_${ASSEMBLY}/${REGION}_bridge_${WINDOW}_chrom_name.vg
		  ./vgRenamePath.sh ${REGION}_${ASSEMBLY}/${REGION}_bridge_${WINDOW}_chrom_name.vg  ${CHROM} ref > ${REGION}_${ASSEMBLY}/${REGION}_bridge_${WINDOW}.vg

		  # to side graph sql (-s option because paths dont cover snps)
		  vg2sg ${REGION}_${ASSEMBLY}/${REGION}_bridge_${WINDOW}.vg ${REGION}_${ASSEMBLY}/database_bridge_${WINDOW}.fa ${REGION}_${ASSEMBLY}/database_bridge_${WINDOW}.sql -s
	 done

done

