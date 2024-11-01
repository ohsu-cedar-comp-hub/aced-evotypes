This workflow documentation is for running samples from an s3 style bucket to scratch space on a node (/mnt/scratch) and publishing back to the s3 bucket (bucket -> scratch -> bucket)

# WGS Nextflow Workflows

## Alignment

script: `workflows/alignment.nf`

Realign files using dockstore-cgpmap with Bwa-mem2. Metrics (Alignment Summary, Wgs, and Insert Size) are generated using Picard. 

## Variant

script: `workflows/variant.nf`

Requires tumor/normal files that have been realigned by CGPMAP using Bwa-mem2. Multiple tools are used in this workflow to perform variant calling and allele counting. This includes: 
    
- ds-CGPWGS (caveman, ascat, brass, pindel, and vagrent)
- delly
- gridss
- manta
- strelka
- mutect2
- allele counter for battenberg

# Reference Files and md5sums

All of the reference files used in the alignment and variant workflows are provided here. See the note below for info regarding downloading these files.

- `core_ref_GRCh38_hla_decoy_ebv.tar.gz`                `6448a15bcc8f91271b1870a3ecfcf630`
- `bwa_idx_GRCh38_hla_decoy_ebv_bwamem2.tar.gz`         `5e1652381e73285d17bc348998fe7612`
- `VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz`  `876657ce8d4a6dd69342a9467ef8aa76`
- `SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz`  `03ac504f1a2c0dbe34ac359a0f8ef690`
- `CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz`      `90a100d06dbde243c6e7e11e6a764374`
- `qcGenotype_GRCh38_hla_decoy_ebv.tar.gz`              `1956e28c1ff99fc877ff61e359e1020c`
- `core_ref_GRCh38_hla_decoy_ebv/genome.fa`             `4bcba916b619324b132d2626c923b5ac`
- `human.hg38.excl.tsv`                                 `55f6e74bf725c1bba89b3e7908cb915f`
- `af-only-gnomad.hg38.vcf.gz`                          `a4209be7fb4b5a5a8d3b778132cb7401`
- `hg38.even.intervals`                                 `0648d246e717be9abf0239cf0c8154c7`
- `small_exac_common_3.hg38.vcf.gz`                     `4c75c1755a45c64e8af7784db7fde009`


## Note

Reference files can be downloaded here: https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/

# Containers

## CGPMAP

**Source Code:** https://github.com/cancerit/dockstore-cgpmap

**Docker Image:** https://quay.io/repository/wtsicgp/dockstore-cgpmap?tab=info

**Tool Version:** `3.3.0`

**Docker Pull Command:** `docker pull quay.io/wtsicgp/dockstore-cgpmap:3.3.0`

**SIF:** `dockstore-cgpmap_3.3.0.sif`

## Picard

**Source Code:** https://github.com/broadinstitute/picard

**Picard Jar:** https://github.com/broadinstitute/picard/releases/tag/3.1.1

**Tool Versions:** `Picard 3.1.1` `Java 17`


## CGPWGS

**Source Code:** https://github.com/cancerit/cgpCaVEManWrapper

**Dockerfile Source Code:** https://github.com/cancerit/dockstore-cgpwgs/tree/2.1.0 

**Docker Image:** https://quay.io/repository/wtsicgp/dockstore-cgpwgs?tab=info

**Docker Pull Command:** `docker pull quay.io/wtsicgp/dockstore-cgpwgs`

**Docker Tag:** `2.1.0`

**SIF:** `dockstore-cgpwgs_2.1.0.sif`

## Delly

**Source Code:** https://github.com/dellytools/delly

**Docker Image:** https://hub.docker.com/r/dellytools/delly/tags

**Docker Pull Command:** `docker pull dellytools/delly:v1.2.6`

**Docker Tag:** `v1.2.6`

**SIF:**  `delly_v1.2.6.sif`


## Strelka and Manta

**Source Code:** https://github.com/cancerit/strelka2-manta

**Docker Image:** https://quay.io/repository/wtsicgp/strelka2-manta

**Docker Pull Command:** `docker pull quay.io/wtsicgp/strelka2-manta`

**Docker Tag:** `latest`

**SIF:** `strelka2-manta_latest.sif`


## GRIDSS

**Source Code:** https://github.com/PapenfussLab/gridss

**Docker Image:** https://quay.io/repository/biocontainers/gridss?tab=info

**Docker Pull Command:** `docker pull quay.io/biocontainers/gridss`

**Docker Tag:** `2.13.2--h50ea8bc_3`

**SIF:** `gridss_2.13.2--h50ea8bc_3.sif`

### Note
The Rscript in the sif has a bug. 

Download https://github.com/PapenfussLab/gridss/releases/download/v2.13.2/gridss-2.13.2.tar.gz and alter libgridss.R as indicated (https://github.com/PapenfussLab/gridss/issues/635):

Locate the fie libgridss.R in your conda env folder (e.g. ~/.conda/envs/<my_env>/share/gridss-2.13.2-2)

Replace line 780 in this file
  
Original:

    isbp = str_detect(VariantAnnotation::fixed(vcf)$ALT, "[\\]\\[]")
    
New:

    isbp = str_detect(as.character(VariantAnnotation::fixed(vcf)$ALT), "[\\]\\[]")  
  
Then rerun gridss_somatic_filter.

This assumes that the ALT fields contain a single allele per line, which seems to be the case in my GRIDSS output VCF files.
<br> 
### Setting Up R Environment

To run `gridss_somatic_filter` the following dependencies are needed:

R 4.3.1 with these libraries:

- argparser
- tidyverse
- stringdist
- testthat
- stringr
- StructuralVariantAnnotation
- rtracklayer
- BSgenome package for your reference genome

## GATK

**Documentation:** https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2

**Docker Image:** https://hub.docker.com/r/broadinstitute/gatk

**Docker Pull Command:** `docker pull broadinstitute/gatk:4.5.0.0`

**Docker Tag:** `4.5.0.0`

**SIF:** `gatk_4.5.0.0.sif`
