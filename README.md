# DRAMS
A tool to Detect and Re-Align Mixed-up Samples based on multi-omics data.

## Highlights
* Need at least three omics types to run this tool.
* Sex information is not necessary, but it’s better to have this information.
* Users can set omics priority depending on which omics type they trust more.

## Software requirement
* Linux operating system
* Python 3.x
* PLINK 1.9
* GCTA

## Get started
To run the sample ID realignment procedure:
```bash
python3 run_DRAMS.py --pair=data/genotypes.merge.highlyrelatedpairs.txt --prior=data/omics_priority --nsex=data/samplelist.nominalSex --gsex=data/samplelist.snpSex --output=data/res
```

## Results
For each sample, DRAMS assign a new ID if the sample has been detected as mixed-up. For the results, a table indicating the original and new sample ID for each sample will be generated. Here is the title and meaning of each column:
1. OmicsType (Omics type)
1. OriginalID (Original sample ID)
1. TrueID (New ID assigned by DRAMS)
1. SwitchedOrNot (Did the ID switched or not?)
1. GeneticSex (Genetic sex of the sample)
1. NominalSex_NewID (Nominal sex of the new ID)
1. SexMatchedOrNot (Did the genetic sex and nominal sex of new ID matched or not?)
<div style='display: none'>
In addition to this table, DRAMS generates another table that is formated as input for Cytoscape visualization. Please see the guidance below on how to visualize sample relationships using Cytoscape.
</div>

## Data preparation
### Call genotypes
We recomment to use the same pipeline to call genotypes for sequencing data of any types, such as DNA sequencing, RNA sequencing, ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing). Here are recommended scripts to call genotypes:
```bash
bash scripts/call_genotypes.step1.sh exampleID  # Map reads to reference genome and Base Quality Score Recalibration (BQSR)
bash scripts/call_genotypes.step2.sh  # Call genotypes by GATK HaplotypeCaller
```

### Check sample contamination
We recommend to check sample contamination and remove contaminated samples before run DRAMS. There are two options to check sample contamination. VerifyBamID (https://genome.sph.umich.edu/wiki/VerifyBamID) can be used based on BAM files for sequencing data. Another option is quite straightforward. We provided an AWK script to check sample contamination based on heterozygous rate. We recommend to remove samples with heterozygous rate largely deviated from other samples.
```bash
# Option1 (VerifyBamID):
ls exampleID*.vcf|while read file; do 
    verifyBamID --vcf $file --bam exampleID.bam --out $file.check --verbose --ignoreRG
done
# Option2 (Heterozygous rate):
ls exampleID*.vcf|while read file; do 
    awk 'BEGIN{FS="\t";OFS="\t"}$1!~/^#/{split($9,a,":");for(i in a){if(a[i]=="GQ") GQi=i;if(a[i]=="DP") DPi=i};split($10,a,":");if(a[GQi]<10||a[DPi]<3||a[DPi]>60) next;if($10~/^0\/0/){hom0++}else if($10~/^0\/1/){het++}else if($10~/^1\/1/){hom1++}}END{print "'$file'",hom0,het,hom1,het/(hom0+het+hom1)}' $file
done > heterozygous_rate.txt  # Calculate heterozygous rate for variants with GQ>=10 and 3<=DP<=60.
```

### Infer genetic sex
We recommend to use Plink to infer genetic sex based on X chromosome heterozygosity and Y chromosome call rate.
```bash
plink --vcf exampleID.vcf --make-bed --out exampleID  # Convert VCF file to PLINK file (PLINK 1.9)
plink --bfile exampleID --split-x hg19 --make-bed --out exampleID  # remove X chromosome pseudo-autosomal region
plink --bfile exampleID --check-sex ycount 0.2 0.8 --out exampleID  # sexcheck (Use default 0.2/0.8 F-statistic thresholds temporarily. As the PLINK website noted, the threshold should be determined by eyeballing the distribution of F-estimates)
```

### Estimate genetic relatedness scores
Genetic relatedness scores among all samples in all omics types were estimated by GCTA.
Be noted that sample ID should be formated like this: "OmicsType|SampleID".
```bash
plink --bfile exampleID1 --bmerge exampleID2.bed exampleID2.bim exampleID2.fam --out exampleID.merge  # Merge input files. This step may be run several times if you have multiple input PLINK files.
gcta64 --bfile exampleID --autosome --maf 0.01 --make-grm --out exampleID  # Estimate genetic relatedness by GCTA
```

### Extract highly related sample pairs
The genetic relatedness scores were in bimodal distribution. We provided a script to extract highly related sample pairs based on the distribution. The results include a table of all highly related sample pairs, a density plot of the sample relatedness scores, and a summary table for the highly related sample pairs for each omics type.
```bash
python3 scripts/extract_highly_related_pairs.py --input_prefix=exampleID --output_prefix=exampleID --threshold=0.65 --min_loci=400
```

<div style='display: none'>
### Guidance on visualize highly related sample pairs through Cytoscape
This is not a necessary step for sample ID realignment, but we recommend to visualize the sample relationships to have a better understanding on how the samples switched.
</div>

