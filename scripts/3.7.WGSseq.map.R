# This script is to converts uBAM to Map ready BAM file
# GatherBqsrReports
# # Merge the recalibration reports resulting from by-interval recalibration

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
####### Step 7 #########
########################

gatk_path = "/home/boss_lab/mgupta/GATK/gatk-4.4.0.0/gatk"
mem_size_gb = 4
command_mem_gb = ceiling(mem_size_gb) - 1

# Set up parallel backend with desired number of cores
# Adjust the value of 'ncores' to specify the number of parallel processes
ncores <- 32
cl <- makeCluster(ncores)
registerDoParallel(cl)

# GatherBqsrReports = function(i){
for (i in 2:nrow(files)){
  print(files$sample)
  input_bqsr_reports <- c()
  for (j in 1:3366) {
    temp_input <- paste0(files$dir[i], files$sample[i], j, ".recal_data.csv")
    input_bqsr_reports <- unique(c(input_bqsr_reports, temp_input))
  }
  output_report_filename = paste0(files$dir[i],files$sample[i],".combined.recal_data.csv")
  command <- paste0(gatk_path,
                    " --java-options", " '-Xms", command_mem_gb, "G'",
                    " GatherBQSRReports",
                    " -I ", paste(input_bqsr_reports, collapse = " -I "),
                    " -O ", output_report_filename)
  
  
  # Execute the command using system()
  # print(command)
  system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
  # return(command)
}


# Iterate over i in parallel
results = foreach(i = 2:nrow(files), .combine = c) %dopar% {
  GatherBqsrReports(i)
}


