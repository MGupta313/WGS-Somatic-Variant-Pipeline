# This script is to converts uBAM to Map ready BAM file
# SamToFastqAndBwaMem
# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment

# Read in bistools
source("/home/boss_lab/Apps/bitbucket/bistools/ESB_bisTools.R") # Magic

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
####### Step 1 #########
########################

# Map reads to reference
refPath = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/pipeline/ref/"

# Set the R variables needed for the command-line
ref_fasta <- paste0(refPath,"Homo_sapiens_assembly38.fasta")
compression_level <- 5
mem_size_gb = 14
command_mem_gb = ceiling(mem_size_gb/2)
gotc_path <- "/home/boss_lab/Apps/anaconda3/share/picard-2.23.8-0/" #picard path

for (i in 1:nrow(files)) {
  
  time_start = Sys.time()
  
  print(paste(files$sample[i], time_start))

  input_bam <- paste0(files$dir[i],files$uBAM[i])
  bwa_path <- "/usr/bin/bwa"
  bwa_commandline <- paste("bwa mem -K 100000000 -p -v 3 -t 16 -Y", ref_fasta)
  output_bam_basename <- paste0(files$dir[i],files$sample[i],".unmerged")

  # Build the command
  command <- paste0("java -Dsamjdk.compression_level=", compression_level, " -Xms", command_mem_gb, "G -jar ", gotc_path, "picard.jar",
                    " SamToFastq INPUT=", input_bam, " FASTQ=/dev/stdout INTERLEAVE=true NON_PF=true | ",
                    bwa_commandline, " /dev/stdin - 2> >(tee ", output_bam_basename, ".bwa.stderr.log >&2) | ",
                    "samtools view -1 - > ", output_bam_basename, ".bam")
  
  # Execute the command using system()
  # print(command)
  system(paste0("/bin/bash -c '", command, "'"))
  
  print(paste("Made unmerged.bam for", files$sample[i]))
}



