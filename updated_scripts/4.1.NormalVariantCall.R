# Step 1. Run Mutect2 in tumor-only mode for each normal sample.
# Call variants in the normal sample

projectDir = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/run2/" # Replace, must have '/' at the end

# Set working directory
homeDir = paste0(projectDir, "analysis/1.PON/")
setwd(homeDir)

# Manifest file
filesDir = paste0(projectDir, "pipeline/") # Replace, if necessary
filesFile = "WGSseq.sample.manifest.txt" # Replace, if necessary
files = read.table(paste0(filesDir, filesFile), sep = "\t", header = T, as.is = T)

normal = files$sample[files$group2 == "Dump"]
normal_dir = files$dir[files$group2 == "Dump"]
tumor = files$sample[files$group2 == "Ro.pos"]
tumor_dir = files$dir[files$group2 == "Ro.pos"]

normal_bam = paste0(normal,".bam")
tumor_bam = paste0(tumor,".bam")

gatk_path = "/home/boss_lab/mgupta/GATK/gatk-4.4.0.0/gatk"
ref_fasta = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/pipeline/ref/Homo_sapiens_assembly38.fasta"

for (i in 1:length(normal_bam)){
  sample = normal[i]
  input = paste0(normal_dir[i], normal_bam[i])
  output = paste0(projectDir,"analysis/1.PON/",normal[i],"_pon.vcf.gz")
  command <- paste0(gatk_path, 
                    " Mutect2 -R ", ref_fasta,
                    " -I ", input,
                    " -max-mnp-distance 0 -O ", output)
  print(command)
  # system(command)
}

# Step 2. Create a GenomicsDB from the normal Mutect2 calls.
# Combined all normal vcfs into a db

cmd2 = paste(gatk_path, "GenomicsDBImport --genomicsdb-workspace-path pon_db -R", ref_fasta, "-V", output, "-V /home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/analysis/PanelOfNormals/SSRR001_pon.vcf.gz -V /home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/analysis/PanelOfNormals/SSRR004_pon.vcf.gz -V /home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/analysis/PanelOfNormals/SSRR007_pon.vcf.gz -V /home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/analysis/PanelOfNormals/SSRR010_pon.vcf.gz -L /home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/analysis/PanelOfNormals/wgs_calling_regions.hg38.interval_list")

system(cmd2)

# Step 3. Combine the normal calls using CreateSomaticPanelOfNormals.
command_mem = 7
gnomad = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/pipeline/ref/af-only-gnomad.hg38.vcf.gz"

cmd3 = paste0(gatk_path, " --java-options '-Xmx7g' CreateSomaticPanelOfNormals -R ", ref_fasta, " --germline-resource ", gnomad, " -V gendb:///home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/run2/analysis/1.PON/pon_db -O /home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/run2/analysis/1.PON/pon.vcf.gz")

system(cmd3, ignore.stdout = TRUE, ignore.stderr = TRUE)

