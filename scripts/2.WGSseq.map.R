# This script is to convert fastq to uBAM

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

for (i in 1:nrow(files)) {
  
  time_start = Sys.time()
  
  print(paste(files$sample[i], time_start))
  # fq1file = paste0(files$dir[i],files$fqMate1[i])
  # fq1 <- readLines(fq1file, n = 1)
  
  gatk = "/home/boss_lab/mgupta/GATK/gatk-4.4.0.0/gatk"
  fastq_1 = paste0(files$dir[i],files$fqMate1[i])
  fastq_2 = paste0(files$dir[i],files$fqMate2[i])
  output = paste0(files$dir[i],files$sample[i], ".unmapped.bam")
  readgroup_name = files$sample[i]
  sample_name = files$sample[i]
  library_name = paste0("Illumina-",readgroup_name)
  platform_unit = "H5TTLDSX7"
  run_date = "2023-05-26"
  platform = "Illumina"
  sequencing_center = "Novogene"
  
  print("Converting Fastq to uBAM")
  # print(paste(gatk, "--java-options '-Xmx6g' FastqToSam --FASTQ",
  #               fastq_1, "--FASTQ2", fastq_2, "--OUTPUT", output,
  #               "--READ_GROUP_NAME", readgroup_name, "--SAMPLE_NAME", sample_name,
  #               "--LIBRARY_NAME",library_name, "--PLATFORM_UNIT",platform_unit, "--RUN_DATE",
  #               run_date,"--PLATFORM Illumina --SEQUENCING_CENTER Novogene"))
  
  system(paste(gatk, "--java-options '-Xmx6g' FastqToSam --FASTQ",
               fastq_1, "--FASTQ2", fastq_2, "--OUTPUT", output,
               "--READ_GROUP_NAME", readgroup_name, "--SAMPLE_NAME", sample_name,
               "--LIBRARY_NAME",library_name, "--PLATFORM_UNIT",platform_unit, "--RUN_DATE",
               run_date,"--PLATFORM Illumina --SEQUENCING_CENTER Novogene"))
  
  files$uBAM[i] = paste0(files$sample[i], ".unmapped.bam")
  files$timeTaken[i] = difftime(Sys.time(), time_start, units = "hours")
  print(paste("Time taken for sample", files$sample[i], "(hrs):", files$timeTaken[i]))
  
  # Update manifest file and note that this sample finished mapping by setting 'include' to FALSE
  files$include[i] = FALSE
  write.table(files, file = paste0(filesDir, filesFile), sep = "\t", row.names = F, quote = F)

}
  
print("All conversion done")
  
  
  
