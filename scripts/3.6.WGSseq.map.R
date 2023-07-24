# This script is to converts uBAM to Map ready BAM file
# BaseRecalibrator
# Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel

# Read in bistools
# source("/home/boss_lab/Apps/bitbucket/bistools/ESB_bisTools.R") # Magic
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
####### Step 6 #########
########################

gatk_path = "/home/boss_lab/mgupta/GATK/gatk-4.4.0.0/gatk"
refPath = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/pipeline/ref/"
ref_fasta <- paste0(refPath,"Homo_sapiens_assembly38.fasta")
mem_size_gb = 6
command_mem_gb = ceiling(mem_size_gb) - 2
dbSNP_vcf = paste0(refPath,"Homo_sapiens_assembly38.dbsnp138.vcf")
known_indels_sites_VCFs = c("/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/pipeline/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz", "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/pipeline/ref/Homo_sapiens_assembly38.known_indels.vcf.gz")

sequence_grouping = read.table("sequence_grouping.txt", sep = "\t", header = F, as.is = T)

# Set up parallel backend with desired number of cores
# Adjust the value of 'ncores' to specify the number of parallel processes
ncores <- 32
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Define a function to be executed in parallel
baseRecalibrator = function(i,j){
  sequence_group_interval = sequence_grouping$V1[j] #subgroup
  
  input_bam = paste0(files$dir[i],files$sample[i],".aligned.duplicate_marked.sorted.bam")
  recalibration_report_filename = paste0(files$dir[i],files$sample[i],j, ".recal_data.csv")
  
  # Build the command string
  command <- paste0(gatk_path,
                    " --java-options", " '-Xms", command_mem_gb, "G'",
                    " BaseRecalibrator",
                    " -R ", ref_fasta,
                    " -I ", input_bam,
                    " --use-original-qualities",
                    " -O ", recalibration_report_filename,
                    " --known-sites ", dbSNP_vcf,
                    " --known-sites ", paste(known_indels_sites_VCFs, collapse = " --known-sites "),
                    " -L ", paste(sequence_group_interval, collapse = " -L "))
  
  # Execute the command using system()
  system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
  return(command)
}

# Iterate over i in parallel
# For every sample, run over each sequence subgroup
results <- foreach (i = 1:nrow(files), .combine = c) %:%
  foreach (j = 1:length(sequence_grouping$V1), .combine = c) %dopar% {
    baseRecalibrator(i, j)
  }




