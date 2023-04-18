IN_BAM=$1
OUT_STATS=$2
samtools stats ${IN_BAM} > ${OUT_STATS}
