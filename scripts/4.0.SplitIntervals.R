# Create a Mutect2 panel of normals
# Interval lists define subsets of genomic regions, sometimes even just individual positions in the genome.
# This tool takes in intervals via the standard arguments of IntervalArgumentCollection and 
# splits them into interval files for scattering. 
# The resulting files contain equal number of bases.

gatk_path = "/home/boss_lab/mgupta/GATK/gatk-4.4.0.0/gatk"
small_task_mem = 4
command_mem= small_task_mem * 1000 - 500
ref_fasta = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/pipeline/ref/Homo_sapiens_assembly38.fasta"
intervals = "/home/boss_lab/Projects/Boss_Projects/WGS.Somatic.SLE/analysis/PanelOfNormals/wgs_calling_regions.hg38.interval_list"
scatter_count = 24
min_contig_size = "1000000"
split_intervals_extra_args = paste0(" --dont-mix-contigs --min-contig-size ",min_contig_size)

system("mkdir interval-files")

command = paste0(gatk_path, " --java-options '-Xmx",command_mem,"m' SplitIntervals",
                 " -R ", ref_fasta,
                 " -L ", intervals,
                 " -scatter ", scatter_count,
                 " -O interval-files",
                 split_intervals_extra_args)
system(command)
#system("cp interval-files/*.interval_list .")

