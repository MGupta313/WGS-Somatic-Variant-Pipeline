# This script is to converts uBAM to Map ready BAM file
# MarkDuplicates
# Mark duplicate reads to avoid counting non-independent observations

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
####### Step 3 #########
########################

# Setting constants
gatk_path = "/home/boss_lab/mgupta/GATK/gatk-4.4.0.0/gatk"
compression_level <- 5
mem_size_gb = 7.5
command_mem_gb = ceiling(mem_size_gb) - 2

# Set up parallel backend with desired number of cores
# Adjust the value of 'ncores' to specify the number of parallel processes
ncores <- 32
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Define a function to be executed in parallel
markDups = function(i){
  time_start = Sys.time()
  
  print(paste(files$sample[i], time_start))
  base_file_name = paste0(files$dir[i],files$sample[i])
  input_bams = paste0(files$dir[i],files$sample[i],".aligned.unsorted.bam")
  output_bam_basename = paste0(files$dir[i],files$sample[i],".aligned.unsorted.duplicates_marked")
  metrics_filename = paste0(files$dir[i],files$sample[i],".duplicate_metrics")
  
  # Build the command string
  command <- paste0(gatk_path,
                    " --java-options '-Dsamjdk.compression_level=", compression_level, " -Xms", command_mem_gb, "G'",
                    " MarkDuplicates",
                    " --INPUT ", paste(input_bams, collapse = " --INPUT "),
                    " --OUTPUT ", output_bam_basename, ".bam",
                    " --METRICS_FILE ", metrics_filename,
                    " --VALIDATION_STRINGENCY SILENT",
                    " --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500",
                    " --ASSUME_SORT_ORDER 'queryname'",
                    " --CREATE_MD5_FILE true")
  
  # Execute the command using system()
  system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
  return(command)
}

# Iterate over i in parallel
results = foreach(i = 1:nrow(files), .combine = c) %dopar% {
  if (i != 3) {
    time_start = Sys.time()
    print(paste(files$sample[i], time_start))
    markDups(i)
  }
}

# Clean up parallel backend
# stopCluster(cl)
# registerDoSEQ()  # Switch back to sequential execution if needed


# for (i in 1:nrow(files)) {
#   
#   time_start = Sys.time()
#   
#   print(paste(files$sample[i], time_start))
#   base_file_name = paste0(files$dir[i],files$sample[i])
#   input_bams = paste0(files$dir[i],files$sample[i],".aligned.unsorted.bam")
#   output_bam_basename = paste0(files$dir[i],files$sample[i],".aligned.unsorted.duplicates_marked")
#   metrics_filename = paste0(files$dir[i],files$sample[i],".duplicate_metrics")
# 
#   # Build the command string
#   command <- paste0(gatk_path,
#                     " --java-options '-Dsamjdk.compression_level=", compression_level, " -Xms", command_mem_gb, "G'",
#                     " MarkDuplicates",
#                     " --INPUT ", paste(input_bams, collapse = " --INPUT "),
#                     " --OUTPUT ", output_bam_basename, ".bam",
#                     " --METRICS_FILE ", metrics_filename,
#                     " --VALIDATION_STRINGENCY SILENT",
#                     " --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500",
#                     " --ASSUME_SORT_ORDER 'queryname'",
#                     " --CREATE_MD5_FILE true")
#   
#   # Execute the command using system()
#   system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
#   
# }
