from gwf import Workflow
import os, sys
import math
from glob import glob

project_name = "CoastSeq"

gwf = Workflow(defaults={"account": "edna"}) 

libraries = [x for x in glob("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/tmp/*") if os.path.isdir(x)]

#Summing sequence tables of all libraries
input_files = []
for library_root in libraries:
    library_id = os.path.basename(library_root)   
    input_files.append("tmp/{}/seqtab_RDS".format(library_id))
    input_files.append("tmp/{}/seqtab.nochim_RDS".format(library_id))
 
output_files = []
    
output_files.append("results/seqtab_Both")
output_files.append("results/seqtab.nochim_Both")
output_files.append("results/DADA2_raw.table")
output_files.append("results/DADA2_raw.otus")
output_files.append("results/DADA2_nochim.table")
output_files.append("results/DADA2_nochim.otus")
    
gwf.target(
   name="sum_libraries_{}".format(project_name), # f'sum_libraries_{project_name}' is equivalent
   inputs=input_files,
   outputs=output_files,
   cores=32,
   memory="128g",
   walltime="4:00:00",
 ) << """
   Rscript ./scripts/sum_libraries.r tmp/ results/
   """.format() 
