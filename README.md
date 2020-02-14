# DRAMS
A tool to Detect and Re-Align Mixed-up Samples for Integrative Studies of Multi-omics Data.

## Highlights
* Need at least three types of omics data to run this tool.
* We used a logistic regression model followed by a modified topological sorting algorithm to correct sample mix-ups based on data relationships of multi-omics. User can also train there own logistic regression model.
* Sex information is not necessary, but having both genetics-based and reported sexes will help identify true IDs.
* Users can set omics priority to indicate their confidence for the correctness of each omics type.

## Software requirement
* Linux operating system
* Python 3.x
* PLINK 1.9
* GCTA

## Get started
Users can use the following script to realign sample IDs, with the input of highly related data pairs, genetics-based and reported sexes, and omics priority. This script firstly estimate possible switch directions for each mismatched data pair based on a logistic regression model. Then, it connects all the highly related data pairs and produces multiple independent networks. A topological sorting algorithm is applied to identify the highest confidence IDs for data IDs in each network. Users can use "--stringent" to run DRAMS in stringent mode. In the stringent mode, we discarded the data groups with less than three data or with no shared IDs (i.e. all data IDs in a group are different), as they are almost unlikely to be corrected.

```bash
python3 run_DRAMS.py --pair=data/genotypes.merge.highlyrelatedpairs.txt --prior=data/omics_priority --nsex=data/samplelist.reportedSex --gsex=data/samplelist.snpSex --output=data/res
```

## Input file examples
1. **Highly related data pairs** (Tab-separated). Please see [document below](#Estimate-genetic-relatedness-scores) for detailed procedure to create this file.

| OmicsType1 | SampleID1 | OmicsType2 | SampleID2 | Relatedness | Match |
| ---------- | --------- | ---------- | --------- | ----------- | ----- |
| omicsA | S1 | omicsB | S1 | 0.977 | Y |
| omicsA | S2 | omicsB | S2 | 0.972 | Y |
| omicsA | S3 | omicsC | S2 | 0.985 | N |

2. **Omics type priority file** (Tab-separated). An integer started from 1 in the second column to denote the user’s confidence for the correctness of each omics type (1 means the highest confidence). The omics types listed in one line separated by comma will be considered as one omics type in the logistic regression.

| OmicsType | Priority |
| ---------- | --------- |
| omicsA | 1 |
| omicsB | 2 |
| omicsC,omicsD | 3 |

3. **Reported sex file** (Tab-separated).

| SampleID | ReportedSex |
| ---------- | --------- |
| S1 | M |
| S2 | F |

4. **Genetics-based sex file** (Tab-separated). Please see [document below](#Infer-genetics-based-sexes) for recommended guideline to estimate genetics-based sexes.

| SampleID | GeneticSex |
| ---------- | --------- |
| S1 | M |
| S2 | F |

## Output file
For each data, DRAMS assign a new ID if the data has been detected as a mix-up. For the results, a table indicating the original and new sample ID for each data will be created. Here is the title and meaning of each column:
1. OmicsType (Omics type)
1. OriginalID (Original sample ID)
1. TrueID (New ID assigned by DRAMS)
1. SwitchedOrNot (Sample ID changed or not)
1. GeneticSex (Genetics-based sex of the data)
1. ReportedSex_NewID (Reported sex of the new ID)
1. SexMatchedOrNot (Genetics-based sex and reported sex matched or not after correcting IDs)

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

### Infer genetics-based sexes
For DNA-based data (such as WGS, ATAC-Seq), we recommend using Plink to infer genetics-based sexes based on X chromosome heterozygosity and Y chromosome call rate. As the PLINK website noted, the threshold should be determined by eyeballing the distribution of F-estimates). Here is the example code of Plink:
```bash
plink --vcf exampleID.vcf --make-bed --out exampleID  # Convert VCF file to PLINK file (PLINK 1.9)
plink --bfile exampleID --split-x hg19 --make-bed --out exampleID  # remove X chromosome pseudo-autosomal region
plink --bfile exampleID --check-sex ycount 0.2 0.8 --out exampleID  # sexcheck (Use default 0.2/0.8 F-statistic thresholds temporarily. As the PLINK website noted, the threshold should be determined by eyeballing the distribution of F-estimates)
```
For RNA-Seq data, we recommend to infer genetics-based sexes according to XIST gene expression levels. In our study, we taked samples with XIST expression larger than 2 (TPM, Transcripts Per Kilobase Million) as females.

### Estimate genetic relatedness scores
Genetic relatedness scores among data of different omics types were estimated by GCTA.
Please be noted that sample ID should be formated like this: "OmicsType|SampleID".
```bash
plink --bfile exampleID1 --merge-list mergelist.txt --make-bed --out exampleID.merge  # Merge input files, each file per line for the merge list.
gcta64 --bfile exampleID --autosome --maf 0.01 --make-grm --out exampleID  # Estimate genetic relatedness by GCTA
```

### Extract highly related data pairs
Basically, the genetic relatedness scores were in bimodal distribution. The threshold to distinguish highly related data pairs from random pairs should be determined by eyeballing the distribution. Based on our BrainGVEX data, we set the threshold as 0.65 by default. We provided a script to extract highly related data pairs from GCTA results. The results include a table of all highly related sample pairs, histograms of the genetic relatedness scores, a Cytoscape input file for visualizing sample relationships, and a summary table for the highly related data pairs for each omics type.
```bash
python3 scripts/extract_highly_related_pairs.py --input_prefix=exampleID --output_prefix=exampleID --threshold=0.65  --plot
```

## Guideline of visualizing data relationships using Cytoscape
The Cytoscape input file named "exampleID.highlyrelatedpairs.cytoscape.txt" has been generated in previous step. Users can visualize the data relationships among all the omics types using Cytoscape. Here is the steps in detail:
1. Load the Cytoscape software (recommended version > 3.0).
1. Click File -> Import -> Network -> File.
1. Select the file "exampleID.highlyrelatedpairs.cytoscape.txt" in the file chooser dialog.
1. Define the first column as Source node and the second column as Target node. Click OK.
1. We recommend to use different edge types to indicate matched and mismatched data pairs. From the Control Panel on the left, select Style tag -> Edge -> Line Type -> Mapping. Choose the "interaction" column, then select Parallel Lines for "Match" and Solid for "Mismatch".

Here is several examples of data relationships visualization:
![GitHub Logo](/images/SampleRelationExample.png)

## Train logistic regression model
We have trained the logistic regression model using our BrainGVEX data. The training results were used as default parameters. We also provide an function for users to train there model using other data if they have. The script was shown below:
```bash
python3 run_DRAMS.py --pair=training/trainingSet.highlyrelatedpairs --prior=training/trainingSet.omicsPriority --nsex=training/trainingSet.reportedSex --gsex=training/trainingSet.snpSex --train=training/trainingSet.highConfDirections --output=training/trainingSet
```

