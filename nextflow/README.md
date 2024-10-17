This workflow documentation is for running samples from an s3 style bucket to scratch space on a node (/mnt/scratch) and publishing back to the s3 bucket (bucket -> scratch -> bucket)

# WGS Nextflow Workflows

alignment.nf
    - realign files by running dockstore-cgpmap with Bwa-mem2
    - generate alignment metrics using picard-3.1.0

variant.nf
    - runs ds-CGPWGS (caveman, ascat, brass, pindel, and vagrent), delly, gridss, manta, strelka, mutect2, battenberg

## Reference Files and md5sums
- `core_ref_GRCh38_hla_decoy_ebv.tar.gz`                `6448a15bcc8f91271b1870a3ecfcf630`
- `bwa_idx_GRCh38_hla_decoy_ebv_bwamem2.tar.gz`         `5e1652381e73285d17bc348998fe7612`
- `VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz`  `876657ce8d4a6dd69342a9467ef8aa76`
- `SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz`  `03ac504f1a2c0dbe34ac359a0f8ef690`
- `CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz`      `90a100d06dbde243c6e7e11e6a764374`
- `qcGenotype_GRCh38_hla_decoy_ebv.tar.gz`              `1956e28c1ff99fc877ff61e359e1020c`
- `core_ref_GRCh38_hla_decoy_ebv/genome.fa`
- `human.hg38.excl.tsv`                                 `55f6e74bf725c1bba89b3e7908cb915f`
- `BSgenome.Hsapiens.UCSC.hg38`
- `af-only-gnomad.hg38.vcf.gz`                          `a4209be7fb4b5a5a8d3b778132cb7401`
- `hg38.even.intervals`                                 `0648d246e717be9abf0239cf0c8154c7`
- `small_exac_common_3.hg38.vcf.gz`                     `4c75c1755a45c64e8af7784db7fde009`


## Note

Instructions on how to download reference files for GRCh38 can be found here: https://github.com/cancerit/dockstore-cgpmap/wiki/1.-Reference-files#grch38

Reference files can also be downloaded here: https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/
