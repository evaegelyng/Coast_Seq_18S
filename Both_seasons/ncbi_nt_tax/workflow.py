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
output_files.append("results/{}_SWARM3nc_stats".format(project_name))
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

###Split fasta file into K parts
#NB! Should add removal of old index file, as this seems to not be overwritten
def splitter(inputFile, K=99):
    inputs = [inputFile]
    outputs = ["tmp/split/split.log.txt"]
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '1:00:00'
    }
    spec = '''
    eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate metabar_2021
    seqkit split -O tmp/split/ {inputFile} -p {K} -2
    echo "hello" > tmp/split/split.log.txt
    '''.format(inputFile=inputFile, K=K)
    return inputs, outputs, options, spec

#####blast a single k-th file
def blaster(k, outFolder):
    inputFasta = 'tmp/split/{}_SWARM_seeds.part_'.format(project_name)+'{:0>3d}'.format(k)+'.fasta'
    inputs = [inputFasta]
    outBlast = outFolder + '/blast.' + str(k) + '.blasthits'
    outLog = outFolder + '/blast.' + str(k) + '.txt'
    outputs = [
      outBlast,
      outLog
    ]
    options = {
        'cores': 2,
        'memory': '32g',
        'walltime': '4:00:00'
    }
    spec = '''
    eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate metabar_2021
    export BLASTDB=/faststorage/project/eDNA/blastdb/nt_220927/
    mkdir -p {out}
    echo "RUNNING THREAD {k} BLAST"
    blastn -db /faststorage/project/eDNA/blastdb/nt_220927/nt -negative_gilist /faststorage/project/eDNA/KR21/OBB/18S/Ole_18S_Negative_GI_list.gi -max_target_seqs 500 -num_threads 4 -outfmt "6 std qlen qcovs staxid" -out {outBlast} -qcov_hsp_perc 90 -perc_identity 90 -query {inputFasta}
    echo "hello" > {outLog}
    echo "hello" > {outLog}
    echo "DONE THREAD {k}"
    '''.format(out=outFolder, k=k, inputFasta=inputFasta, outBlast=outBlast, outLog=outLog)
    return inputs, outputs, options, spec

def taxonomy(taxonomyFolder, blastFolder, k):
    inputFile = blastFolder + '/blast.' + str(k) + '.blasthits'
    inputs = [inputFile , blastFolder + '/blast.' + str(k) + '.txt']
    summaryFile = taxonomyFolder + '/summary.' + str(k) + '.txt'
    outputFile = taxonomyFolder + '/taxonomy.' + str(k) + '.txt'
    outputs = [summaryFile, outputFile]
    options = {
        'cores': 1,
        'memory': '32g',
        'walltime': '4:00:00'
    }
    
    spec = '''
    eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate metabar_2021
    mkdir -p {taxonomyFolder}
    # Check if blast file is empty
    if [ `cat {inputFile} | wc -l` != 0 ]
    then
      Rscript scripts/taxonomy_v0_1_1_EES.R {inputFile} {summaryFile} {outputFile}
    else
      touch {outputFile}
    fi
    if grep -q "Query coverage is less than 100% for all hits" ".gwf/logs/taxonomy_" + str(k) + ".stdout"
    then
      touch {outputFile}
    if grep -q "Sequence identity is less than 90% for all hits" ".gwf/logs/taxonomy_" + str(k) + ".stdout"
    then
      touch {outputFile}
    '''.format(taxonomyFolder=taxonomyFolder, inputFile=inputFile, summaryFile=summaryFile, outputFile=outputFile) 
    return inputs, outputs, options, spec

inputName = "results/{}_SWARM_seeds.fasta".format(project_name)

gwf.target_from_template( 'split', splitter(inputFile=inputName) )

parts=glob('tmp/split/{}_SWARM_seeds.part*.fasta'.format(project_name))
K=len(parts)
                                                                
for k in range(1,K+1):
  gwf.target_from_template( 'blaster_{}'.format(k), blaster(k=k, outFolder='tmp/blast') )
  gwf.target_from_template( 'taxonomy_{}'.format(k), taxonomy(taxonomyFolder='tmp/taxonomy', blastFolder='tmp/blast', k=k) )

### Combine all the small taxonomical summary files into one large file

input_files = glob('tmp/taxonomy/summary*.txt')
 
output_file = 'results/{}_summary.txt'.format(project_name)
    
gwf.target(
   name="combine_summary_{}".format(project_name),
   inputs=input_files,
   outputs=output_file,
   cores=1,
   memory="1g",
   walltime="00:10:00",
 ) << """
    head -n1 tmp/taxonomy/summary.1.txt > results/{project_name}_summary.txt
    for fname in tmp/taxonomy/summary*.txt
    do
        tail -n +2 $fname >> results/{project_name}_summary.txt
    done
   """.format(project_name=project_name)    
      
### Combine all the small taxonomical classfication files into one large file

input_files = glob('tmp/taxonomy/taxonomy*.txt')
 
output_file = 'results/{}_classified.tsv'.format(project_name)
    
gwf.target(
   name="combine_taxonomy_{}".format(project_name),
   inputs=input_files,
   outputs=output_file,
   cores=1,
   memory="1g",
   walltime="00:10:00",
 ) << """
    head -n1 tmp/taxonomy/taxonomy.1.txt > results/{project_name}_classified.tsv
    for fname in tmp/taxonomy/taxonomy*.txt
    do
        tail -n +2 $fname >> results/{project_name}_classified.tsv
    done
   """.format(project_name=project_name)         
