# This script is to converts uBAM to Map ready BAM file
# GatherBamFiles
# Merge the recalibrated BAM files resulting from by-interval recalibration

# Read in bistools
source("/home/boss_lab/Apps/bitbucket/bistools/ESB_bisTools.R") # Magic
library(foreach)
library(doParallel)

projectDir = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/" # Replace, must have '/' at the end

# Set working directory
homeDir = paste0(projectDir, "pipeline/")
setwd(homeDir)

# Manifest file
filesDir = homeDir # Replace, if necessary
filesFile = "WGSseq.sample.manifest.txt" # Replace, if necessary
files = read.table(paste0(filesDir, filesFile), sep = "\t", header = T, as.is = T)

# Set output directory
outDir = paste0(projectDir, "data/")

########################
####### Step 9 #########
########################

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs

gatk_path = "/home/boss_lab/mgupta/GATK/gatk-4.4.0.0/gatk"
mem_size_gb = 3
command_mem_gb = ceiling(mem_size_gb) - 1
compression_level = 5

# Set up parallel backend with desired number of cores
# Adjust the value of 'ncores' to specify the number of parallel processes
ncores <- 32
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Define a function to be executed in parallel
GatherBamFiles = function(i){
  
  input_bams <- c()
  for (j in 1:3366) {
    temp_input <- paste0(files$dir[i],files$sample[i],".",j,".aligned.duplicates_marked.recalibrated.bam")
    input_bams <- unique(c(input_bams, temp_input))
  }
  output_bam_basename = paste0(files$dir[i],files$sample[i])
  
  # Build the command string
  command <- paste0(gatk_path, 
                    " --java-options '-Dsamjdk ", compression_level, " -Xms", command_mem_gb, "G'",
                    " GatherBamFiles ",
                    "--INPUT ", paste(input_bams, collapse = " --INPUT "), " ",
                    "--OUTPUT ", output_bam_basename, ".bam ",
                    "--CREATE_INDEX true ",
                    "--CREATE_MD5_FILE true")
  
  # Execute the command using system()
  # system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
  return(command)
}

# Iterate over i in parallel
results = foreach(i = 2:nrow(files), .combine = c) %dopar% {
  GatherBamFiles(i)
}






