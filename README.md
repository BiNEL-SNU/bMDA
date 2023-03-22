bMDA
## Requirements
- required software : BWA-MEM, SAMtools, Picard Toolkit, Genome Analysis Toolkit (GATK v3.7-0), BEDTools, VarScan2, MuTect (v1.1.4), MuTect2, DELLY, Manta, GRIDSS2, SURVIVOR, AnnotSV, Circos
- dependent R package : DNACopy, aCGH, copynumber, ape, ggtree, ChromoMap (v4.1.1), KataegisPortal
- Operating system : Ubuntu
- installing all required software would finish within 2 hours.

## Installation guide
```shell
git clone https://github.com/BiNEL-SNU/bMDA.git
```

## Demo samples
- demo data can be found in example
- Processing demo sample will complete within 30min

## Instructions for use
### SNV calling
```shell
bash SNV_call/SNV_call.sh
```
Output: matrix designating detected SNVs

### SV calling
```shell
bash SV_call/SV_call.sh
Rscript SV_call/analysis of SV_190422.R
```
Output: matrix designating detected SVs


### Kataegis calling
```shell
Rscript Kataegis_call/kataegis_event_190422.R
```
Output: Rainfall plot in pdf format
