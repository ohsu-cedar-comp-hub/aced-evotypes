# Use an official R base image with R version 4.3.1
FROM rocker/r-ver:4.3.1

# Install system dependencies required for R package installations
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libgit2-dev \
    gfortran \
    wget bzip2 \
    && rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# Install Miniconda (or Conda, if necessary) for managing bioconda packages
RUN wget -qO- https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh > miniconda.sh && \
    bash miniconda.sh -b -f -p /root/miniconda3 && \
    rm miniconda.sh && \
    /root/miniconda3/bin/conda clean --all -y

# Set up environment variables for Conda
ENV PATH=/root/miniconda3/bin:$PATH

# Install a specific Python version if necessary
RUN /root/miniconda3/bin/conda install -y python=3.10

# Install cancerit-allelecount from Bioconda with the specified version
RUN /root/miniconda3/bin/conda install -y -c bioconda -c conda-forge cancerit-allelecount=4.3.0

# Clean up Conda cache to reduce the image size
RUN /root/miniconda3/bin/conda clean --all -y

# Install specific R packages and versions
RUN R -e "install.packages('doParallel', repos='https://cran.r-project.org', version='1.0.17')" \
    && R -e "install.packages('iterators', repos='https://cran.r-project.org', version='1.0.14')" \
    && R -e "install.packages('foreach', repos='https://cran.r-project.org', version='1.5.2')" \
    && R -e "install.packages('optparse', repos='https://cran.r-project.org', version='1.7.4')" \
    && R -e "install.packages('getopt', repos='https://cran.r-project.org', version='1.20.4')" \
    && R -e "install.packages('codetools', repos='https://cran.r-project.org', version='0.2.19')"

# Ensure the container starts with R (optional)
CMD ["R"]
