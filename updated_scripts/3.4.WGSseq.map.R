# This script is to converts uBAM to Map ready BAM file
# SortAndFixTags
# Sort BAM file by coordinate order and fix tag values for NM and UQ

# Read in bistools
source("/home/boss_lab/Apps/bitbucket/bistools/ESB_bisTools.R") # Magic
library(foreach)
library(doParallel)

projectDir = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/run2/" # Replace, must have '/' at the end

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
####### Step 4 #########
########################

gatk_path = "/home/boss_lab/mgupta/GATK/gatk-4.4.0.0/gatk"
refPath = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/pipeline/ref/"
ref_fasta <- paste0(refPath,"Homo_sapiens_assembly38.fasta")
compression_level <- 5
mem_size_gb = 10
command_mem_gb_sort = ceiling(mem_size_gb) - 1
command_mem_gb_fix = ceiling((mem_size_gb - 1)/10)

bwa_version = "0.7.17-r1188"
bwa_commandline <- paste("bwa mem -K 100000000 -p -v 3 -t 16 -Y", ref_fasta)

# Set up parallel backend with desired number of cores
# Adjust the value of 'ncores' to specify the number of parallel processes
ncores <- 8
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Sort BAM file by coordinate order and fix tag values for NM and UQ

# Define a function to be executed in parallel
sortAndFixTags = function(i){
  time_start = Sys.time()
  
  input_bam = paste0(files$dir[i],files$sample[i],".aligned.unsorted.duplicates_marked.bam")
  output_bam_basename = paste0(files$dir[i],files$sample[i],".aligned.duplicate_marked.sorted")
  
  # Build the first command string
  command1 <- paste0(gatk_path,
                     " --java-options '-Dsamjdk.compression_level=", compression_level, " -Xms", command_mem_gb_sort, "G'",
                     " SortSam",
                     " --INPUT ", input_bam,
                     " --OUTPUT /dev/stdout",
                     " --SORT_ORDER 'coordinate'",
                     " --CREATE_INDEX false",
                     " --CREATE_MD5_FILE false")

  # Build the second command string
  command2 <- paste0(gatk_path,
                     " --java-options '-Dsamjdk.compression_level=", compression_level, " -Xms", command_mem_gb_fix, "G'",
                     " SetNmMdAndUqTags",
                     " --INPUT /dev/stdin",
                     " --OUTPUT ", output_bam_basename, ".bam",
                     " --CREATE_INDEX true",
                     " --CREATE_MD5_FILE true",
                     " --REFERENCE_SEQUENCE ", ref_fasta)

  # Combine the two commands using a pipe
  full_command <- paste(command1, "|", command2)
  
  # Execute the command using system()
  system(full_command, ignore.stdout = TRUE, ignore.stderr = TRUE)
  return(full_command)
}

# Iterate over i in parallel
results = foreach(i = 1:nrow(files), .combine = c) %dopar% {
  sortAndFixTags(i)
}

# Clean up parallel backend
stopCluster(cl)
registerDoSEQ()  # Switch back to sequential execution if needed


  

