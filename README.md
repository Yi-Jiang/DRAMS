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

## Input files
1. Highly related data pairs (Tab-separated). Please see [document below](#Estimate-genetic-relatedness-scores) for detailed procedure to create this file.
    OmicsType1	SampleID1	OmicsType2	SampleID2	Relatedness	Match
    omics2	S1	omics1	S1	0.9778299331665039	Y
    omics2	S2	omics1	S2	0.9724648594856262	Y
    omics2	S3	omics1	S3	1.0035074949264526	Y

## Results
For each sample, DRAMS assign a new ID if the sample has been detected as mixed-up. For the results, a table indicating the original and new sample ID for each sample will be generated. Here is the title and meaning of each column:
1. OmicsType (Omics type)
1. OriginalID (Original sample ID)
1. TrueID (New ID assigned by DRAMS)
1. SwitchedOrNot (Did the ID switched or not?)
1. GeneticSex (Genetic sex of the sample)
1. NominalSex_NewID (Nominal sex of the new ID)
1. SexMatchedOrNot (Did the genetic sex and nominal sex of new ID matched or not?)

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
plink --bfile exampleID1 --merge-list mergelist.txt --make-bed --out exampleID.merge  # Merge input files, each file per line for the merge list.
gcta64 --bfile exampleID --autosome --maf 0.01 --make-grm --out exampleID  # Estimate genetic relatedness by GCTA
```

### Extract highly related sample pairs
The genetic relatedness scores were in bimodal distribution. We provided a script to extract highly related sample pairs based on the distribution. The results include a table of all highly related sample pairs, histograms of the sample relatedness scores, a Cytoscape input file for visualizing sample relationships, and a summary table for the highly related sample pairs for each omics type.
```bash
python3 scripts/extract_highly_related_pairs.py --input_prefix=exampleID --output_prefix=exampleID --threshold=0.65  --plot
```

### Guideline of visualizing sample relationships using Cytoscape
The Cytoscape input file named "exampleID.highlyrelatedpairs.cytoscape.txt" has been generated in previous step. Users can visualize the sample relationships among all the -omics types using Cytoscape. Here is the steps in detail:
1. Load the Cytoscape software (recommended version > 3.0).
1. Click File -> Import -> Network -> File.
1. Select the file "exampleID.highlyrelatedpairs.cytoscape.txt" in the file chooser dialog.
1. Define the first column as Source node and the second column as Target node. Click OK.
1. We recommend to use different edge types to indicate matched and mismatched sample pairs. From the Control Panel on the left, select Style tag -> Edge -> Line Type -> Mapping. Choose the "interaction" column, then select Parallel Lines for "Match" and Solid for "Mismatch".

Here is several examples of sample relationships visualization:
![GitHub Logo](/images/SampleRelationExample.png)


