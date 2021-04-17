#!/bin/bash
echo "trimming"
# To improve the quality of fastq reads, different commands are used to perform trimming stratigies. Trimming is important to help in allignment with the reference genome and remove all bad bases and any contaminating vector or adapter sequences from the reads as well as filter out poor quality score reads.
# We have to run our raw reads through FastQC to assess the quality of our sequencing reads (R1 and R2 reads). Now we are going to improve the quality of our reads, by trimming off any "bad" bases using the following trimmomatic commands (for paired-end fastq files):

trimmomatic PE \
  -threads 4 \
  -phred33 \
  $1 $2 \
  -baseout ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data \
  ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10
  TRAILING:25 MINLEN:50

echo "fastQC"
# This section is for assessing the quality of our data using fastQC tools:

fastqc -t 4 /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_1P \
	/home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_2P

echo "Alignment"
#Alignment with bwa and bwa mem:
# Fllowing trimming, the fllowing steps are used to align sequence reads to the reference genome using alignment method BWA.
# Run the commands: $ bwa and $ bwa mem for read allignment
# We need to index the genome with `bwa index' command.
# For read group (RG) information, run bwa mem with RG information using the following commands:

bwa mem -t 4 -v 1 -R '@RG\tID:11V6WR1.111.HD1375ACXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera\tDT:2017-02-28\tPU:11V6WR1' -I 250,50  ~/ngs_course/dnaseq/data/reference/hg19.fa.gz ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_1P ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_2P > ~/ngs_course/dnaseq/data/aligned_data/NGS0001.sam

echo "Duplicate_marking"
# The commands below are used for post alignment quality control (QC) and filtering. There are 2 steps: 
(1) Using Picard tools to mark duplicated reads and (2) Filtering bam data
# We use Picard tools to mark duplicated reads. This tool examines aligned records in the sam and bam file, where can locate duplicate molecules.
# We have to mark duplicates first and then filtering the bam file data as below.

picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt

echo "post_alignment_filtering"
# Filtering bam data dependent on mapping quality and bitwise flags using samtools:
# Filtering is based on different criteria:
#Minimum MAPQ quality score : 20 -Filter on bitwise flag: yes a. Skip alignments with any of these flag bits.

samtools view -F 1796 -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam
samtools index NGS0001_sorted_filtered.bam

echo "variant_calling"
# This section for variants calling using freebayes software.
# FreeBayes is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms),
# This process requires two inputs: a FASTA reference sequence, and a bam-format alignment file by reference position.
# The output reports positions which it finds polymorphic in variant calling file (VCF) format. 
# To convert to text format the reference (as required by samtools faidx) and index it with samtools faidx, we have to use the following commands.
# The following commands used for calling variants with Freebayes tools and then compress the output (the results are variants in vcf file):

freebayes --bam ~/ngs_course/dnaseq/data/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference ~/ngs_course/dnaseq/data/reference/hg19.fa --vcf ~/ngs_course/dnaseq/results/NGS0001.vcf

echo "variant_filtering"
# This can be performed via this command: $ vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1".
# That means: "QUAL > 1: removes horrible sites QUAL / AO > 10 : additional contribution of each obs should be 10 log units 
# The command for filtering variants is here:

vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ~/ngs_course/dnaseq/results/NGS0001.vcf.gz > ~/ngs_course/dnaseq/results/NGS0001_filtered.vcf

echo "variant_annotation_annovar"
# To run the ANNOVAR table function, we generate the following command and the results and output in csv format: 

./table_annovar.pl ~/ngs_course/dnaseq/results/NGS0001_filtered_chr22.avinput humandb/ -buildver hg19 -out ~/ngs_course/dnaseq/results/NGS0001_filtered_chr22 -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout

echo "variant_annotation_snpEFF"
# To generate the variants calling format file using snpEFF, we used this commands:
java -Xmx8g -jar /home/ubuntu/snpEFF/snpEff/snpEff.jar hg19 ~/ngs_course/dnaseq/results/NGS0001_filtered_chr22.vcf.gz > 
/home/ubuntu/ngs_course/dnaseq/results/NGS0001_filtered.ann.vcf
