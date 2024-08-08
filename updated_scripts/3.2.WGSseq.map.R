# This script is to converts uBAM to Map ready BAM file
# MergeBamAlignment
# Merge original input uBAM file with BWA-aligned BAM file

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
####### Step 2 #########
########################

# Merge original input uBAM file with BWA-aligned BAM file

# Set the R variable needed for the command-line
gatk_path = "/home/boss_lab/mgupta/GATK/gatk-4.4.0.0/gatk"
refPath = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/pipeline/ref/"
ref_fasta <- paste0(refPath,"Homo_sapiens_assembly38.fasta")
compression_level <- 5
mem_size_gb = 4
command_mem_gb = ceiling(mem_size_gb) - 1
bwa_version = "0.7.17-r1188"
bwa_commandline <- paste("bwa mem -K 100000000 -p -v 3 -t 16 -Y", ref_fasta)

bash_ref_fasta <- ref_fasta

# Set up parallel backend with desired number of cores
# Adjust the value of 'ncores' to specify the number of parallel processes
ncores <- 32
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Define a function to be executed in parallel
mergeFunc = function(i){
  
  time_start = Sys.time()
  print(paste(files$sample[i], time_start))
  
  aligned_bam = paste0(files$dir[i],files$sample[i],".unmerged.bam")
  unmapped_bam = paste0(files$dir[i],files$uBAM[i])
  output_bam_basename = paste0(files$dir[i],files$sample[i],".aligned.unsorted")
  
  
  # Build the command string
  command <- paste0(gatk_path," --java-options '-Dsamjdk.compression_level=", compression_level, " -Xms", command_mem_gb, "G\'",
                    " MergeBamAlignment",
                    " --VALIDATION_STRINGENCY SILENT",
                    " --EXPECTED_ORIENTATIONS FR",
                    " --ATTRIBUTES_TO_RETAIN X0",
                    " --ALIGNED_BAM ", aligned_bam,
                    " --UNMAPPED_BAM ", unmapped_bam,
                    " --OUTPUT ", output_bam_basename, ".bam",
                    " --REFERENCE_SEQUENCE ", ref_fasta,
                    " --PAIRED_RUN true",
                    " --SORT_ORDER 'unsorted'",
                    " --IS_BISULFITE_SEQUENCE false",
                    " --ALIGNED_READS_ONLY false",
                    " --CLIP_ADAPTERS false",
                    " --MAX_RECORDS_IN_RAM 2000000",
                    " --ADD_MATE_CIGAR true",
                    " --MAX_INSERTIONS_OR_DELETIONS -1",
                    " --PRIMARY_ALIGNMENT_STRATEGY MostDistant",
                    " --PROGRAM_RECORD_ID 'bwamem'",
                    " --PROGRAM_GROUP_VERSION '", bwa_version, "'",
                    " --PROGRAM_GROUP_COMMAND_LINE '", bwa_commandline, "'",
                    " --PROGRAM_GROUP_NAME 'bwamem'",
                    " --UNMAPPED_READ_STRATEGY COPY_TO_TAG",
                    " --ALIGNER_PROPER_PAIR_FLAGS true",
                    " --UNMAP_CONTAMINANT_READS true")
  
  # Execute the command using system()
  system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
  return(command)
}

# Iterate over i in parallel
results = foreach(i = 1:nrow(files), .combine = c) %dopar% {
  mergeFunc(i)
}
