inputBAM="-I example1.realigned.recal.bam -I example2.realigned.recal.bam -I example3.realigned.recal.bam"

# call genotypes
java -Xmx20g -jar GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar -l INFO -R refgenome.hg19.fa -T HaplotypeCaller -nct 6 --dbsnp All_20180423.vcf.gz -o vcf/sampleall.raw.vcf -L $interval $inputBAM
