## List all input BAM files here
inputBAM="-I example1.realigned.recal.bam -I example2.realigned.recal.bam -I example3.realigned.recal.bam"

## call genotypes
# GATK can be downloaded from https://software.broadinstitute.org/
# All_20180423.vcf.gz downloaded from ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/
java -Xmx20g -jar GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar -l INFO -R refgenome.hg19.fa -T HaplotypeCaller -nct 5 --dbsnp All_20180423.vcf.gz -o vcf/sampleall.raw.vcf $inputBAM

## For large data, recommended to call genotypes for each small region separately.
#java -Xmx20g -jar GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar -l INFO -R refgenome.hg19.fa -T HaplotypeCaller -nct 5 --dbsnp All_20180423.vcf.gz -o vcf/sampleall.raw.vcf -L chr1:1-5000000 $inputBAM
