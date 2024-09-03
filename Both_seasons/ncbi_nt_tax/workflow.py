from gwf import Workflow
import os, sys
import math
from glob import glob

project_name = "COSQ"

gwf = Workflow(defaults={"account": "edna"}) 

#Using Mjolnir pipeline to perform OTU clustering.
# Run add_count_fa.r before running this workflow

input_files = []

input_files.append("../results/no_sing_cleaned_otu_table_ASV_wise.txt")
input_files.append("../results/COSQ.fa")

output_files = []
 
output_files.append("results/{}_SWARM_seeds.fasta".format(project_name))
output_files.append("results/{}_SWARM13nc_stats".format(project_name))
output_files.append("results/{}_SWARM_output".format(project_name))
output_files.append("results/{}_non_singleton_motu_list.txt".format(project_name))
output_files.append("results/{}_SWARM_output_counts.tsv".format(project_name))

gwf.target(
            name="odin_{}".format(project_name),
            inputs=input_files,
            outputs=output_files,
            cores=16,
            memory="164g",
            walltime="4:00:00",            
        ) << """
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate dnoise_2
            cd results
            Rscript ../scripts/odin_220224.r ~/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/ {project_name}
        """.format(project_name=project_name) 
