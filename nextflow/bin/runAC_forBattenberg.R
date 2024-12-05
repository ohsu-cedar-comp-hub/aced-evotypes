# load in libraries
library("doParallel")
library("iterators")
library("foreach")
library("optparse")
library("compiler")
library("getopt")
library("codetools")
library("parallel")
library("tools")

# get options
option_list = list(
  make_option(c("-t", "--tumourname"), type="character", default=NULL, help="Samplename of the tumour", metavar="character"),
  make_option(c("-n", "--normalname"), type="character", default=NULL, help="Samplename of the normal", metavar="character"),
  make_option(c("--tb"), type="character", default=NULL, help="Tumour BAM file", metavar="character"),
  make_option(c("--nb"), type="character", default=NULL, help="Normal BAM file", metavar="character"),
  make_option(c("--cpu"), type="numeric", default=8, help="The number of CPU cores to be used by the pipeline (Default: 16)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Directory where output will be written", metavar="character"),
  # make_option(c("-w", "--workdir"), type="character", default=NULL, help="Working directory path", metavar="character"),
  make_option(c("-p", "--allelesprefix"), type="character", default=NULL, help="g1000allelesprefix", metavar="character")
  )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# verify options
print("t:")
print(opt$tumourname)
print("n:")
print(opt$normalname)
print("tb:")
print(opt$tb)
print("nb:")
print(opt$nb)
print("cpu:")
print(opt$cpu)
print("o:")
print(opt$output)
print("p:")
print(opt$allelesprefix)

tumourname = opt$tumourname
normalname = opt$normalname
normalbam = opt$nb
tumourbam = opt$tb
RUN_DIR = opt$output
g1000allelesprefix = opt$allelesprefix

print("")
print(paste("RUN_DIR is: ", RUN_DIR))
print("")

# number threads is set to 16 but this can be changed
nthreads = opt$cpu

# set working directory to run directory
setwd(RUN_DIR)

# want to count bases with quality thresholds
min_base_qual = 20
min_map_qual = 35


chrom_names <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10","11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X")

# Initialize parallel backend
cl <- makeCluster(nthreads)
registerDoParallel(cl)


#' Obtain allele counts for 1000 Genomes loci through external program alleleCounter -v 4.3.0
#'
#' @param bam.file A BAM alignment file on which the counter should be run.
#' @param output.file The file where output should go.
#' @param g1000.loci A file with 1000 Genomes SNP loci.
#' @param min.base.qual The minimum base quality required for it to be counted (optional, default=20).
#' @param min.map.qual The minimum mapping quality required for it to be counted (optional, default=35).
#' @param allelecounter.exe A pointer to where the alleleCounter executable can be found (optional, default points to $PATH).
getAlleleCounts = function(bam.file, output.file, g1000.loci, min.base.qual=20, min.map.qual=35, allelecounter.exe="alleleCounter") {
  cmd = paste(allelecounter.exe,
              "-b", bam.file,
              "-l", g1000.loci,
              "-o", output.file,
              "-m", min.base.qual,
              "-q", min.map.qual,
	   "--dense-snps")
    system(cmd, wait=T)
	}

alleleCounting = function(chrom_names, tumourbam, normalbam, tumourname, normalname, g1000allelesprefix, min_base_qual, min_map_qual,   nthreads) {
    # Export the getAlleleCounts function to the parallel workers
    foreach::foreach(i=1:length(chrom_names), .export=c("getAlleleCounts")) %dopar% {
        getAlleleCounts(bam.file=tumourbam,
                        output.file=paste(tumourname,"_alleleFrequencies_chr", chrom_names[i], ".txt", sep=""),
                        g1000.loci=paste(g1000allelesprefix, chrom_names[i], ".txt", sep=""),
                        min.base.qual=min_base_qual,
                        min.map.qual=min_map_qual)
      
        getAlleleCounts(bam.file=normalbam,
                        output.file=paste(normalname,"_alleleFrequencies_chr", chrom_names[i], ".txt",  sep=""),
                        g1000.loci=paste(g1000allelesprefix, chrom_names[i], ".txt", sep=""),
                        min.base.qual=min_base_qual,
                        min.map.qual=min_map_qual)
    }
}


# Call the  allelecounting function
alleleCounting(chrom_names, tumourbam, normalbam, tumourname, normalname, g1000allelesprefix, min_base_qual, min_map_qual, nthreads) 

# Stop the cluster
stopCluster(cl)

# want to check that all files have been made and are not too small
# Minimum size in bytes
min_size <- 7494

# Get a list of files with 'alleleFrequencies' in their names
files <- list.files(RUN_DIR, pattern = "alleleFrequencies", full.names = TRUE)

# Check if the number of files is 46
if (length(files) != 46) {
  cat("Job failed: Number of files is not 46.\n")
  quit(status = 1)
}

# Check the size of each file
for (file in files) {
  file_size <- file.info(file)$size
  if (is.na(file_size) || file_size < min_size) {
    cat("Job failed: File", file, "is smaller than", min_size, "bytes.\n")
    quit(status = 1)
  }
}

cat("Job done\n")
# Output package versions
sessionInfo()