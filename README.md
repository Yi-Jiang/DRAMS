# DRAMS
A tool to Detect and Re-Align Mixed-up Samples based on multi-omics data.

## Highlights
* Need at least three omics types to run this tool.
* Sex information is not necessary, but it’s better to have this information.
* Users can set omics priority depending on which omics type they trust more.

## System requirement
DRAMS can be run under linux environment with Python3 installed. No installation required.

## Get started
To run the sample ID realignment procedure:
```bash
python3 run_DRAMS.py -a relatedness.highlyRelate.filter -b omics_priority -c samplelist --coef="0,4.41,8.94,0.19" --prefix=relatedness.highlyRelate
```

## Results


## Data preparation
__* Call genotypes__


__* Check sample contamination__


__* Infer genetic sex__


__* Estimate sample relatedness scores__


__* Extract highly related sample pairs__


__* Summarise results__



