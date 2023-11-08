###############################################################
#Completion of alignment and variant detection using E.coli K12
###############################################################
#Construction of the reference genome sequence of E.coli K12
mkdir E.coli
cd E.coli
mkdir fasta
cd fasta
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gzip -dc GCF_000005845.2_ASM584v2_genomic.fna.gz > E.coli_K12_MG1655.fa
samtools faidx E.coli_K12_MG1655.fa

#Download the sequencing data of E.coli K12
cd ..
mkdir fastq
cd fastq
fastq-dump --split-files SRR1770413

#Quality control
#Quality control of raw sequencing data
mkdir fastqc_out_dir
fastqc *.fastq -o fastqc_out_dir/
#Excise sequencing adapter sequences and low-quality reads
java -jar /mnt/c/Users/yutc/Desktop/wgs_practice/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog logfile SRR1770413_1.fastq SRR1770413_2.fastq out.SRR1770413_1.fastq out.trim.SRR1770413_1.fastq out.SRR1770413_2.fastq out.trim.SRR1770413_2.fastq ILLUMINACLIP:/mnt/c/Users/yutc/Desktop/wgs_practice/bin/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:50

#Alignment Locate sequencing data to the reference genome
#Before formal alignment, the FM-index (alignment index) required for BWA alignment needs to be constructed for the reference sequence.
cd ..
cd fasta
bwa index E.coli_K12_MG1655.fa
#bwa completes comparison
cd ..
bwa mem -t 4 -R '@RG\tID:foo\tPL:illumina\tSM:E.coli_K12' /mnt/c/Users/yutc/Desktop/wgs_practice/input/E.coli/fasta/E.coli_K12_MG1655.fa /mnt/c/Users/yutc/Desktop/wgs_practice/input/E.coli/fastq/SRR1770413_1.fastq /mnt/c/Users/yutc/Desktop/wgs_practice/input/E.coli/fastq/SRR1770413_2.fastq > E_coli_K12.sam
samtools view -bS E_coli_K12.sam >E_coli_K12.bam
samtools sort -@ 4 -m 1G -O bam -o E_coli_K12.sorted.bam E_coli_K12.bam && echo "** BAM sort done"
rm -f E_coli_K12.bam

##label PCR duplication
gatk MarkDuplicates -I E_coli_K12.sorted.bam -O E_coli_K12.sorted.markdup.bam -M E_coli_K12.sorted.markdup_metrics.txt && echo "** markdup done **"
rm -f E_coli_K12.bam
rm -f E_coli_K12.sorted.bam

#Create a comparison index file
samtools index E_coli_K12.sorted.markdup.bam && echo "** index done **"

#Mutation detection
#generate a. dict file for the reference sequence of E.coli K12, which can be completed by calling the CreateSequenceDictionary module
gatk CreateSequenceDictionary -R fasta/E.coli_K12_MG1655.fa -O E.coli_K12_MG1655.dict && echo "** dict done **"
#Using the HaplotypeCaller module of GATK 4.0
#Generate a GVCF for each sample first, and then use GenotypeGVCFs to perform joint calls on these GVCFs
##1. Generate intermediate file gvcf
gatk HaplotypeCaller -R fasta/E.coli_K12_MG1655.fa --emit-ref-confidence GVCF -I E_coli_K12.sorted.markdup.bam -O E_coli_K12.g.vcf && echo "** gvcf done **"
##2. Mutation detection through GVCF
gatk GenotypeGVCFs -R fasta/E.coli_K12_MG1655.fa -V E_coli_K12.g.vcf -O E_coli_K12.vcf && echo "** vcf done **"

#Compress this VCF using bgzip and build an index for it using tabix for future analysis.
#1 compression
bgzip -f E_coli_K12.vcf
#2. Building a tabix index
tabix -p vcf E_coli_K12.vcf.gz

#Variation quality control
# Use SelectVariants to select SNP
gatk SelectVariants -select-type SNP -V E_coli_K12.vcf.gz -O E_coli_K12.snp.vcf.gz
# Perform hard filtering for SNPs
gatk VariantFiltration -V E_coli_K12.snp.vcf.gz --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "PASS" -O E_coli_K12.snp.filter.vcf.gz
# Use SelectVariants to select Indel
gatk SelectVariants -select-type INDEL -V E_coli_K12.vcf.gz -O E_coli_K12.indel.vcf.gz
# Filter for Indel
gatk VariantFiltration -V E_coli_K12.indel.vcf.gz --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "PASS" -O E_coli_K12.indel.filter.vcf.gz
# Re-merge the filtered SNP and Indel
gatk MergeVcfs -I E_coli_K12.snp.filter.vcf.gz -I E_coli_K12.indel.filter.vcf.gz -O E_coli_K12.filter.vcf.gz
# Delete useless intermediate files
rm -f E_coli_K12.snp.vcf.gz* E_coli_K12.snp.filter.vcf.gz* E_coli_K12.indel.vcf.gz* E_coli_K12.indel.filter.vcf.gz*

