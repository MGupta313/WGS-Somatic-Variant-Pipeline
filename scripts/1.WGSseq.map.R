# This script is to convert fastq to unampped bam file format
# Read in bistools
source("/home/boss_lab/Apps/bitbucket/bistools/ESB_bisTools.R") # Magic

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

threads = 32
library_name = "" #"Solexa-NA12878"
platform_unit = "" #"H06HDADXX130110.2.ATCACGAT"
sequencing_center = "Novogene" #"BI"
run_date = "2023-05-24" #give today's date



# paste("./gatk-4.4.0.0/gatk --java-options '-Xmx6g' FastqToSam --FASTQ",
#       fastq1, "--FASTQ2", fastq2, "--OUTPUT" , output, "--READ_GROUP_NAME", read_grp_name,"--SAMPLE_NAME",        sample_name, "--LIBRARY_NAME", library_name, "--PLATFORM_UNIT", platform_unit, "--RUN_DATE", run_date,       "--PLATFORM illumina", "--SEQUENCING_CENTER", sequencing_center
#       )


for (i in 1:nrow(files)) {
  
  if (!files$include[i]) { next }
  
  time_start = Sys.time()
  print(paste(files$sample[i], time_start))
  
  print("Copy files to data dir") # Record where the data was copied from first
  sampleDir = paste0(outDir, files$sample[i], "/")
  if (!file.exists(paste0(sampleDir, files$fqMate1[i])) || !file.exists(paste0(sampleDir, files$fqMate2[i]))) {
    files$fqMate1_input[i] = paste0(files$dir[i], stringr::str_split_1(files$fqMate1[i], "\\|"), collapse = "|")
    files$fqMate2_input[i] = paste0(files$dir[i], stringr::str_split_1(files$fqMate2[i], "\\|"), collapse = "|")
  }
  mvFiles(files$fqMate1[i], files$dir[i], sampleDir)
  files$dir[i] = mvFiles(files$fqMate2[i], files$dir[i], sampleDir)
  
  print("Format and concatenate fastq files")
  files$fqMate1[i] = formatFastq(files$fqMate1[i], files$dir[i], paste0(files$sample[i], "_1"), threads = threads)
  files$fqMate2[i] = formatFastq(files$fqMate2[i], files$dir[i], paste0(files$sample[i], "_2"), threads = threads)
  print("combining of fastq files done")
  
  # print("Converting fastq to bam")
  
  files$timeTaken[i] = difftime(Sys.time(), time_start, units = "hours")
  print(paste("Time taken for sample", files$sample[i], "(hrs):", files$timeTaken[i]))
  
  # Update manifest file and note that this sample finished mapping by setting 'include' to FALSE
  files$include[i] = FALSE
  write.table(files, file = paste0(filesDir, filesFile), sep = "\t", row.names = F, quote = F)
}

# Update manifest file
files$include = TRUE
write.table(files, file = paste0(filesDir, filesFile), sep = "\t", row.names = F, quote = F)
