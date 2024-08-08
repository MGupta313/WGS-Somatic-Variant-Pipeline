# Create somatic panel of normals
gatk_path = "/home/boss_lab/mgupta/GATK/gatk-4.4.0.0/gatk"
gnomad = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/pipeline/ref/af-only-gnomad.hg38.vcf.gz"
ref_fasta = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/pipeline/ref/Homo_sapiens_assembly38.fasta"
pon_db = "gendb:///home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/run2/analysis/1.PON/pon_db"
pon_vcf = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/run2/analysis/1.PON/pon.vcf.gz"

cmd = paste(gatk_path, "--java-options '-Xmx3500g' CreateSomaticPanelOfNormals -R", ref_fasta, "--germline-resource", gnomad, "-V", pon_db, "-O", pon_vcf)

system(cmd)
