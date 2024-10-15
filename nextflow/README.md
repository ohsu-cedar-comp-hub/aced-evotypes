This workflow documentation is for running samples from an s3 style bucket to scratch space on a node (/mnt/scratch) and publishing back to the s3 bucket (bucket -> scratch -> bucket)

# WGS Nextflow Workflows

alignment.nf
    - realign files by running dockstore-cgpmap with Bwa-mem2
    - generate alignment metrics using picard-3.1.0

variant.nf
    - 




## Reference Files and md5sums
- `bwa_idx_GRCh38_hla_decoy_ebv_bwamem2.tar.gz`     `5e1652381e73285d17bc348998fe7612`
- `core_ref_GRCh38_hla_decoy_ebv.tar.gz`            `6448a15bcc8f91271b1870a3ecfcf630`


## Note

Instructions on how to download reference files for GRCh38 can be found here: https://github.com/cancerit/dockstore-cgpmap/wiki/1.-Reference-files#grch38

Reference files can also be downloaded here: https://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/