# Create PON databse
command = paste("gatk GenomicsDBImport -R /home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/pipeline/ref/Homo_sapiens_assembly38.fasta -L /home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/analysis/PanelOfNormals/wgs_calling_regions.hg38.interval_list --genomicsdb-workspace-path pon_db -V /home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/analysis/PanelOfNormals/SSRR001_pon.vcf.gz -V /home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/analysis/PanelOfNormals/SSRR004_pon.vcf.gz -V /home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/analysis/PanelOfNormals/SSRR007_pon.vcf.gz -V /home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/analysis/PanelOfNormals/SSRR010_pon.vcf.gz")

# This script combines the normal calls using CreateSomaticPanelOfNormals

gatk_path = "/home/boss_lab/mgupta/GATK/gatk-4.4.0.0/gatk"

normals_for_pon_vcf_args = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/analysis/PanelOfNormals/pon.args"

#old version:
command <- paste0(gatk_path, 
                  " CreateSomaticPanelOfNormals -V ", normals_for_pon_vcf_args,
                  " -O /home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/analysis/PanelOfNormals/pon.vcf.gz")

system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)

# latest version:
# gatk CreateSomaticPanelOfNormals -R reference.fasta -V gendb://pon_db -O pon.vcf.gz

#github version:
# gatk --java-options "-Xmx~{command_mem}g"  CreateSomaticPanelOfNormals -R ~{ref_fasta} --germline-resource ~{gnomad} \
# -V gendb://pon_db -O ~{output_vcf_name}.vcf ~{create_pon_extra_args}



intervals = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/analysis/PanelOfNormals/wgs_calling_regions.hg38.interval_list"