# This script creates somatic panel of normals
library(foreach)
library(doParallel)

# Set up parallel backend with desired number of cores
# Adjust the value of 'ncores' to specify the number of parallel processes
ncores <- 32
cl <- makeCluster(ncores)
registerDoParallel(cl)

projectDir = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/" # Replace, must have '/' at the end

# Set working directory
homeDir = paste0(projectDir, "pipeline/")
setwd(homeDir)

# Manifest file
filesDir = homeDir # Replace, if necessary
filesFile = "WGSseq.sample.manifest.txt" # Replace, if necessary
files = read.table(paste0(filesDir, filesFile), sep = "\t", header = T, as.is = T)

normal = files$sample[files$group == "Dump_Pos"]
normal_dir = files$dir[files$group == "Dump_Pos"]
tumor = files$sample[files$group == "Ro_Pos_Non_ASC"]
tumor_dir = files$dir[files$group == "Dump_Pos"]

normal_bam = paste0(normal,".bam")
tumor_bam = paste0(tumor,".bam")

gatk_path = "/home/boss_lab/mgupta/GATK/gatk-4.4.0.0/gatk"
ref_fasta = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/pipeline/ref/Homo_sapiens_assembly38.fasta"


CreateSomaticPON = function(i){
# for (i in 1:length(normal_bam)){
  sample = normal[i]
  input = paste0(normal_dir[i], normal_bam[i])
  output = paste0(projectDir,"analysis/PanelOfNormals/",normal[i],"_pon.vcf.gz")
  command <- paste0(gatk_path, 
                    " Mutect2 -R ", ref_fasta,
                    " -I ", input,
                    " -max-mnp-distance 0 -O ", output)
  system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
  return(command)
}

# Iterate over i in parallel
results = foreach(i = 1:length(normal_bam), .combine = c) %dopar% {
  CreateSomaticPON(i)
}


# gatk Mutect2 -R reference.fasta -I normal1.bam -max-mnp-distance 0 -O normal1.vcf.gz
