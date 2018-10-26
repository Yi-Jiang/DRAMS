## Read sample name
sample=$1

## Map reads to reference genome
# BWA can be downloaded from http://bio-bwa.sourceforge.net/
bwa mem -t 5 -R "@RG\tID:foo\tPL:foo\tSM:${sample}" refgenome.hg19.fa ${sample}.r1.fq.gz ${sample}.r2.fq.gz | samtools view -b -S - | samtools sort -m 500000000 -T ${sample}.temp - -o ${sample}.sorted.bam
samtools index ${sample}.sorted.bam

## Split BAM (Recommended for large files)
# samtools view -b -h ${sample}.sorted.bam $chr > ${sample}.sorted.chr$chr.bam
# samtools index ${sample}.sorted.chr$chr.bam

## Mark duplicates 
# Picard-tools can be downloaded from https://sourceforge.net/projects/picard/files/picard-tools/1.119/
java -jar picard-tools-1.119/MarkDuplicates.jar INPUT=${sample}.sorted.bam OUTPUT=${sample}.markdup.bam METRICS_FILE=${sample}.metrics PROGRAM_RECORD_ID=null
samtools index ${sample}.markdup.bam

## Base Quality Score Recalibration (BQSR)
# GATK can be downloaded from https://software.broadinstitute.org/
# All_20180423.vcf.gz downloaded from ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/
java -Xms14g -Xmx18g -jar GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar -T RealignerTargetCreator -R refgenome.hg19.fa -o ${sample}.bam.list -I ${sample}.markdup.bam -nt 5
java -Xms14g -Xmx18g -Djava.io.tmpdir=${sample}.tmpdir -jar GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar -I ${sample}.markdup.bam -R refgenome.hg19.fa -T IndelRealigner -targetIntervals ${sample}.bam.list -o ${sample}.realigned.bam --maxReadsForRealignment 30000 --maxReadsInMemory 1000000
samtools index ${sample}.realigned.bam

java -Xms14g -Xmx18g -jar GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar -l INFO -R refgenome.hg19.fa -T BaseRecalibrator --knownSites:dbsnp,VCF All_20180423.vcf.gz -I ${sample}.realigned.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o ${sample}.recal_data.csv --default_platform illumina -nct 5
java -Xms14g -Xmx18g -jar GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar -l INFO -R refgenome.hg19.fa -I ${sample}.realigned.bam -T PrintReads -o ${sample}.realigned.recal.bam -BQSR ${sample}.recal_data.csv
samtools index ${sample}.realigned.recal.bam
