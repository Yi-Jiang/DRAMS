# DRAMS
A tool to Detect and Re-Align Mixed-up Samples based on multi-omics data.

## Highlights
* Need at least three omics types to run this tool.
* Sex information is not necessary, but it’s better to have this information.
* Users can set omics priority depending on which omics type they trust more.

## System requirement
DRAMS can be run under Linux environment with Python3 installed.

## Get started
To run the sample ID realignment procedure:
```bash
python3 run_DRAMS.py -a relatedness.highlyRelate.filter -b omics_priority -c samplelist --coef="0,4.41,8.94,0.19" --prefix=relatedness.highlyRelate
```

## Results


## Data preparation
### Call genotypes
We recomment to use the same pipeline to call genotypes for sequencing data of any types, such as DNA sequencing, RNA sequencing, ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing). Here are recommended scripts to call genotypes:
```bash
bash Utilities/call_genotypes.step1.sh example  # Map reads to reference genome and Base Quality Score Recalibration (BQSR)
bash Utilities/call_genotypes.step2.sh  # Call genotypes by GATK HaplotypeCaller
```

### Check sample contamination


### Infer genetic sex


### Estimate genetic relatedness scores
Genetic relatedness scores among all samples in all omics types were estimated by GCTA.
```bash
plink --vcf example.vcf --make-bed --out example  # Convert VCF file to PLINK file (PLINK 1.9)
plink --bfile example1 --bmerge example2.bed example2.bim example2.fam --out example.merge  # Merge input files. This step may be run several times if you have multiple input PLINK files.
gcta64 --bfile example --autosome --maf 0.01 --make-grm --out example  # Estimate genetic relatedness by GCTA
```

### Extract highly related sample pairs
The genetic relatedness scores were in bimodal distribution. We provided a script to extract highly related sample pairs based on the distribution. The results include a table of all highly related sample pairs, a density plot of the sample relatedness scores, and a summary table for the highly related sample pairs for each omics type. (Also a formated input for Cytoscape?)
```base
python3 Utilities/extract_highly_related_pairs.py --input_prefix=example --output_prefix=example --threshold=0.65 --min_loci=400
```

### Guidance on visualize highly related sample pairs through Cytoscape
This is not a necessary step for sample ID realignment, but we recommend to visualize the sample relationships to have a better understanding on how the samples switched.


