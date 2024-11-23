# Use the bioconductor/bioconductor_docker:RELEASE_3_17 image as the base. Includes R-4.3.1
FROM bioconductor/bioconductor_docker:RELEASE_3_17

# Install system dependencies required for R package installations
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libgit2-dev \
    gfortran \
    wget bzip2 \
    && rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# Install R packages 
RUN R -e "install.packages('tidyverse', repos='https://cran.r-project.org', version='2.0.0')"
RUN R -e "install.packages('argparser', dependencies=TRUE)"
RUN R -e "install.packages('stringdist', dependencies=TRUE)"
RUN R -e "install.packages('testthat', dependencies=TRUE)"
RUN R -e "install.packages('stringr', dependencies=TRUE)"
RUN R -e "BiocManager::install('StructuralVariantAnnotation', ask=FALSE)"
RUN R -e "BiocManager::install('rtracklayer', ask=FALSE)"
RUN R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38', ask=FALSE)"