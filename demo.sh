
## Preprocess for the four omics data
for i in {1..4}; do
    ## Convert VCF to PLINK
    plink --vcf data/genotypes.omics$i.vcf --make-bed --out data/genotypes.omics$i
    
    ## Add nominal sex information to FAM file
    mv data/genotypes.omics$i.fam data/genotypes.omics$i.fam.bak
    awk 'NR==FNR&&FNR>1{a["omics'$i'|"$1]=$2}NR>FNR{$5=a[$2];print $0}' data/samplelist.nominalSex data/genotypes.omics$i.fam.bak > data/genotypes.omics$i.fam
    
    ## Remove X chromosome pseudo-autosomal region
    plink --bfile data/genotypes.omics$i --split-x hg19 --make-bed --out data/genotypes.omics$i
    
    ## sexcheck (Use default 0.2/0.8 F-statistic thresholds temporarily. As the PLINK website noted, the threshold should be determined by eyeballing the distribution of F-estimates)
    plink --bfile data/genotypes.omics$i --check-sex ycount 0.2 0.8 --out data/genotypes.omics$i
done

## Combine SNP sexes
awk 'BEGIN{OFS="\t";print "OmicsType\tSampleID\tSNP_Sex"}FNR>1{print $2,$4}' data/genotypes.omics*.sexcheck|sed -e's/|/\t/' > data/samplelist.snpSex

## Merge the four omics data
plink --bfile data/genotypes.omics1 --merge-list data/filelist --make-bed --out data/genotypes.merge  # Merge input files

## Estimate genetic relatedness
gcta64 --bfile data/genotypes.merge --autosome --maf 0.01 --make-grm --out data/genotypes.merge  # Estimate genetic relatedness by GCTA

## Extract highly related pairs
python3 scripts/extract_highly_related_pairs.py --input_prefix=data/genotypes.merge --output_prefix=data/genotypes.merge

## Determine true IDs for each samples.
python3 run_DRAMS.py --pair=data/genotypes.merge.highlyrelatedpairs.txt --prior=data/omics_priority --nsex=data/samplelist.nominalSex --gsex=data/samplelist.snpSex --output=data/res
