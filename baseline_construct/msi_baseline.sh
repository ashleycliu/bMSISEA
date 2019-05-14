usage() { echo "Usage: $0 [-r <path to reference genome .fasta/.fa>] [-s <path to baseline construction scripts>] [-o <output path>] [-b <custom_assay_bed>] [-t <output path of msi-h tumor, paired normal coverage file of all covered sites>] [-n <output path of normal/mss coverage file of all covered sites >] [-p <output path of negative plasma coverage file of all covered sites>] [-f <msih_tumor_normal_tumorPercent_file; Header: tumor\tnormal\ttumor_bam\tnormal_bam\tumorPercent>] [-g <mss_normal_list; Header: ID\tbam>] [-i <neg_pla_list; Header: ID\tbam>]" 1>&2; exit 1; }

while getopts ":r:s:o:b:t:n:p:f:g:i:" o; do
    case "${o}" in
        r)
            reference_fasta=${OPTARG}
            ;;
        s)
            script_dir=${OPTARG}
            ;;
        o)
            outdir=${OPTARG}
            ;;
	b)
	    custom_assay_bed=${OPTARG}
            ;;
	t)
	    msihDir=${OPTARG}
            ;;
	n)
            normalDir=${OPTARG}
            ;;
	p)
            negPlaDir=${OPTARG}
            ;;
	f)
	    msih_tumor_normal_tumorPercent_file=${OPTARG}
            ;;
	g)
	    mss_normal_list=${OPTARG}
            ;;
        i)
	    neg_pla_list=${OPTARG}
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${reference_fasta}" ] || [ -z "${outdir}" ] || [ -z "${script_dir}" ] || [ -z "${custom_assay_bed}" ] || [ -z "${msihDir}" ] || [ -z "${normalDir}" ] || [ -z "${negPlaDir}" ] || [ -z "${msih_tumor_normal_tumorPercent_file}" ] || [ -z "${mss_normal_list}" ] || [ -z "${neg_pla_list}" ]; then 
    usage
fi

# 1) Catalog all homo-polymers in your host genome.
${script_dir}/msisensor-raw scan -d ${reference_fasta} -o ${outdir}/genome_microsatellites.list -l 15
paste <(awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$2}' <(sed '1d' ${outdir}/genome_microsatellites.list)) <(awk 'BEGIN{FS="\t";OFS="-"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' <(sed '1d' ${outdir}/genome_microsatellites.list)) > ${outdir}/genome_microsatellites.bed

# 2) Limit the list of microsatellites to those presented in your capture design.
custom_assay_bed_name=`basename ${custom_assay_bed}`
custom_assay_bed_name="${custom_assay_bed_name%.*}"
bedtools intersect -a  ${outdir}/genome_microsatellites.bed -b ${custom_assay_bed} -f 1 -wa > ${outdir}/${custom_assay_bed_name}.microsatellites.bed
cut -f 4 ${outdir}/${custom_assay_bed_name}.microsatellites.bed| awk 'BEGIN{FS="-";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > ${outdir}/${custom_assay_bed_name}.msi.list.tmp
head -1  ${outdir}/genome_microsatellites.list | cat - ${outdir}/${custom_assay_bed_name}.msi.list.tmp > ${outdir}/${custom_assay_bed_name}.msi.list
rm ${outdir}/${custom_assay_bed_name}.msi.list.tmp -f
# 3)  Select microsatellite marker loci and construct baseline files
sed '1d' ${msih_tumor_normal_tumorPercent_file} | while read i; do
    s_tumor=`echo $i|cut -f 1 -d " "`
    s_normal=`echo $i|cut -f 2 -d " "`
    bam_tumor=`echo $i|cut -f 3 -d " "`
    bam_normal=`echo $i|cut -f 4 -d " "`
    ${script_dir}/msisensor-raw msi -d ${outdir}/${custom_assay_bed_name}.msi.list -n ${bam_normal} -t ${bam_tumor} -o ${msihDir}/${s_tumor}_allSite_msi
done 

sed '1d' ${mss_normal_list} | while read i; do
    s_normal=`echo $i|cut -f 1 -d " "`
    bam_normal=`echo $i|cut -f 2 -d " "`
    ${script_dir}/msisensor-raw msi -d ${outdir}/${custom_assay_bed_name}.msi.list -n ${bam_normal} -t ${bam_normal} -o ${normalDir}/${s_normal}_allSite_msi
done 

sed '1d' ${neg_pla_list} | while read i; do
    s_pla=`echo $i|cut -f 1 -d " "`
    bam_pla=`echo $i|cut -f 2 -d " "`
    ${script_dir}/msisensor-raw msi -d ${outdir}/${custom_assay_bed_name}.msi.list -n ${bam_pla} -t ${bam_pla} -o ${negPlaDir}/${s_pla}_allSite_msi
done 

${rscript} ${script_dir}/constructBaseline.R --panelName ${custom_assay_bed_name}  --normalDir ${normalDir} --msihDir ${msihDir} --negplaDir ${negPlaDir} --tumorPercentageFile msih_tumor_normal_tumorPercent_file --targetDepth ${sequencingDepth} --minSupportReads ${minSupportReads}

