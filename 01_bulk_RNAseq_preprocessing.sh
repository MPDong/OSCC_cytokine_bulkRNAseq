#!/bin/bash
# Bulk RNA-seq: QC, trimming, alignment, transcript quantification and counts.
# Run from the project directory; set reference paths to your installation.
base=$(pwd)
sample=$(cat "${base}/name.txt")
raw=${base}/raw
fastqc=${base}/fastqc
fastp=${base}/fastp
fastqc_post=${base}/fastqc_post
ref=${base}/ref
index=${ref}/index_2021/genome_tran
mapped=${base}/mapped
map_sum=${base}/map_sum
mkdir -p "$fastqc" "$fastp" "$fastqc_post" "$mapped" "$map_sum"

# raw/ contains one R1 and one R2 FASTQ per sample, named to match
# ${name}_*_R1_001.fastq.gz and ${name}_*_R2_001.fastq.gz.
# For technical runs of the same library, concatenate matching reads first.
# Example (run separately after setting name and the input directory):
# cat downloads/${name}_*_R1_001.fastq.gz > raw/${name}_merged_R1_001.fastq.gz
# cat downloads/${name}_*_R2_001.fastq.gz > raw/${name}_merged_R2_001.fastq.gz

# FastQC
for name in $sample
do
  fastqc -o "$fastqc" ${raw}/${name}*.fastq.gz
done

# Compile QC reports
multiqc .

# fastp and post-trimming FastQC
for name in $sample
do
  r1=${raw}/${name}_*_R1_001.fastq.gz
  r2=${raw}/${name}_*_R2_001.fastq.gz
  r1_trimmed=${fastp}/${name}_R1_trim.fastq.gz
  r2_trimmed=${fastp}/${name}_R2_trim.fastq.gz

  fastp -i $r1 \
        -I $r2 \
        -o "$r1_trimmed" \
        -O "$r2_trimmed" \
        -j "$fastp/${name}.json" \
        -h "$fastp/${name}.html"

  fastqc -o "$fastqc_post" ${fastp}/${name}*trim.fastq.gz
done

# Alignment and transcript quantification
for name in $sample
do
  hisat2 -p 4 --dta \
    -x "$index" \
    -1 "${fastp}/${name}_R1_trim.fastq.gz" \
    -2 "${fastp}/${name}_R2_trim.fastq.gz" \
    -S "${mapped}/${name}.sam" \
    --summary-file "${map_sum}/${name}_hisat2.log"

  samtools sort -@ 4 -o "${mapped}/${name}.bam" "${mapped}/${name}.sam"
  rm "${mapped}/${name}.sam"

  mkdir -p "${mapped}/${name}"
   stringtie -p 4 -e -B \
	-A ${name}_gene_abund.tab \
	-G ${ref}/Homo_sapiens.GRCh38.102.gtf \
	-o ${mapped}/${name}/${name}.gtf \
	${mapped}/${name}.bam
done

# Assemble the sample/GTF input for prepDE from per-sample GTF directories.
ls "${mapped}"/*/*.gtf | xargs -n 1 basename | cut -d '.' -f 1 > "${mapped}/file1"
ls "${mapped}"/*/*.gtf | xargs -n 1 readlink -f > "${mapped}/file2"
paste "${mapped}/file1" "${mapped}/file2" > "${mapped}/file_path.txt"

# Gene and transcript count matrices
cd "$mapped" || exit 1
python3 "${ref}/prepDE.py3" \
  -i file_path.txt \
  -g GC.csv \
  -t TC.csv
