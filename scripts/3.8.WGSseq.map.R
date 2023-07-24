# This script is to converts uBAM to Map ready BAM file
# ApplyBQSR
# Apply Base Quality Score Recalibration (BQSR) model

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
####### Step 8 #########
########################

# Apply the recalibration model by interval

gatk_path = "/home/boss_lab/mgupta/GATK/gatk-4.4.0.0/gatk"
mem_size_gb = 4
command_mem_gb = ceiling(mem_size_gb) - 1

# Set up parallel backend with desired number of cores
# Adjust the value of 'ncores' to specify the number of parallel processes
ncores <- 32
cl <- makeCluster(ncores)
registerDoParallel(cl)


sequence_grouping = read.table("sequence_grouping.txt", sep = "\t", header = F, as.is = T)
# sequence_group_interval = sequence_grouping$V1 or subgroup # ?????

refPath = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/pipeline/ref/"
ref_dict = paste0(refPath,"Homo_sapiens_assembly38.dict")
ref_fasta = paste0(refPath,"Homo_sapiens_assembly38.fasta")
ref_fasta_index = paste0(refPath,"Homo_sapiens_assembly38.fasta.fai")

# Define a function to be executed in parallel
ApplyBQSR = function(i,j){
  sequence_group_interval = sequence_grouping$V1[j] #subgroup
  
  input_bam = paste0(files$dir[i],files$sample[i],".aligned.duplicate_marked.sorted.bam")
  input_bam_index = paste0(files$dir[i],files$sample[i],".aligned.duplicate_marked.sorted.bai")
  output_bam_basename = paste0(files$dir[i],files$sample[i],".",j,".aligned.duplicates_marked.recalibrated")
  recalibration_report = paste0(files$dir[i],files$sample[i],".combined.recal_data.csv")
  
  # Build the command string
  command <- paste0(gatk_path, 
                    " --java-options", " '-Xms", command_mem_gb, "G'",
                    " ApplyBQSR -R ", ref_fasta, " -I ", input_bam,
                    " -O ", output_bam_basename, ".bam",
                    " -L ", paste(sequence_group_interval, collapse = " -L "),
                    " -bqsr ", recalibration_report,
                    " --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30",
                    " --add-output-sam-program-record",
                    " --create-output-bam-md5",
                    " --use-original-qualities")
  
  # Execute the command using system()
  system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
  return(command)
}

# Iterate over i in parallel
# For every sample, run over each sequence subgroup
results <- foreach (i = 1:nrow(files), .combine = c) %:%
  foreach (j = 1:length(sequence_grouping$V1), .combine = c) %dopar% {
    ApplyBQSR(i, j)
  }















