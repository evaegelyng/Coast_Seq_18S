from gwf import Workflow
import os, sys
import math
from glob import glob

project_name = "CoastSeq"

gwf = Workflow(defaults={"account": "edna"}) 

#Demultiplex

batchfile = "batchfileDADA2.list"

libraries = [x for x in glob("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/base/data/*") if os.path.isdir(x)]

for library_root in libraries:
    library_id = os.path.basename(library_root)
    input_files = glob(library_root + "/*.fq")
    
    with open(os.path.join(library_root, "tags.txt")) as tags_file:
        for line in tags_file:
            output_files = []
            tag_id, fseq, rseq = line.split()

            output_files.append("tmp/{}/DADA2_AS/{}_R1.fastq".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_AS/{}_R2.fastq".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_SS/{}_R1.fastq".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_SS/{}_R2.fastq".format(library_id, tag_id))

            gwf.target(
                name="demultiplex_{}_{}_{}".format(project_name,library_id, tag_id),
                inputs=input_files,
                outputs=output_files,
                cores=1,
                memory="4g",
                walltime="2:00:00",
            ) << """
                mkdir -p tmp/{library_id}
                ./scripts/demultiplex.sh {library_root} tmp/{library_id} {tag_id} {tag_fseq} {tag_rseq} {batchfile}
            """.format(library_root=library_root, library_id=library_id, tag_id=tag_id, tag_fseq=fseq, tag_rseq=rseq, batchfile=batchfile) 
            
#Quality trimming of reads
for library_root in libraries:
    library_id = os.path.basename(library_root)
    
    with open(os.path.join(library_root, "tags.txt")) as tags_file:
        for line in tags_file:
            input_files = []
            tag_id, fseq, rseq = line.split()
    
            input_files.append("tmp/{}/DADA2_AS/{}_R1.fastq".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_AS/{}_R2.fastq".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_SS/{}_R1.fastq".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_SS/{}_R2.fastq".format(library_id, tag_id))
    
            output_files = []
            tag_id, fseq, rseq = line.split()

            output_files.append("tmp/{}/DADA2_AS/filtered/{}_F_filtered.fastq".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_AS/filtered/{}_R_filtered.fastq".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_SS/filtered/{}_F_filtered.fastq".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_SS/filtered/{}_R_filtered.fastq".format(library_id, tag_id))
            
            folderAS="tmp/{}/DADA2_AS/filtered".format(library_id)
            folderSS="tmp/{}/DADA2_SS/filtered".format(library_id)            
            
            gwf.target(
                name="sickle_{}_{}_{}".format(project_name,library_id, tag_id),
                inputs=input_files,
                outputs=output_files,
                cores=1,
                memory="4g",
                walltime="1:00:00",
            ) << """
                mkdir -p {folderAS}
                mkdir -p {folderSS}
                sickle se -l 70 -q 28 -x -t sanger -f {inputASF} -o {outputASF} || echo "Empty file " {outputASF} 
                sickle se -l 70 -q 28 -x -t sanger -f {inputASR} -o {outputASR} || echo "Empty file " {outputASR} 
                sickle se -l 70 -q 28 -x -t sanger -f {inputSSF} -o {outputSSF} || echo "Empty file " {outputSSF} 
                sickle se -l 70 -q 28 -x -t sanger -f {inputSSR} -o {outputSSR} || echo "Empty file " {outputSSR} 
            """.format(folderAS=folderAS,folderSS=folderSS,inputASF=input_files[0],inputASR=input_files[1],inputSSF=input_files[2],inputSSR=input_files[3],outputASF=output_files[0],outputASR=output_files[1],outputSSF=output_files[2],outputSSR=output_files[3])
       
#Matching paired reads
for library_root in libraries:
    library_id = os.path.basename(library_root)
    
    with open(os.path.join(library_root, "tags.txt")) as tags_file:
        for line in tags_file:
            input_files = []
            tag_id, fseq, rseq = line.split()
    
            input_files.append("tmp/{}/DADA2_AS/filtered/{}_F_filtered.fastq".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_AS/filtered/{}_R_filtered.fastq".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_SS/filtered/{}_F_filtered.fastq".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_SS/filtered/{}_R_filtered.fastq".format(library_id, tag_id))
    
            output_files = []
            tag_id, fseq, rseq = line.split()
            
            output_files.append("tmp/{}/DADA2_AS/filtered/matched/{}_F_matched.fastq.gz".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_AS/filtered/matched/{}_R_matched.fastq.gz".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_SS/filtered/matched/{}_F_matched.fastq.gz".format(library_id, tag_id))
            output_files.append("tmp/{}/DADA2_SS/filtered/matched/{}_R_matched.fastq.gz".format(library_id, tag_id))
                        
            folderAS="tmp/{}/DADA2_AS/filtered/matched".format(library_id)
            folderSS="tmp/{}/DADA2_SS/filtered/matched".format(library_id)            
            
            gwf.target(
                name="match_{}_{}_{}".format(project_name,library_id, tag_id),
                inputs=input_files,
                outputs=output_files,
                cores=1,
                memory="8g",
                walltime="2:00:00",
            ) << """
                mkdir -p {folderAS}
                mkdir -p {folderSS}
                ASFfilesize=`stat -c %s {inputASF}`
                if [ $ASFfilesize = 0 ]
                then
                 touch {outputASF} {outputASR}
                fi
                ASRfilesize=`stat -c %s {inputASR}`
                if [ $ASRfilesize = 0 ]
                then
                 touch {outputASF} {outputASR}
                fi
                SSFfilesize=`stat -c %s {inputSSF}`
                if [ $SSFfilesize = 0 ]
                then
                 touch {outputSSF} {outputSSR}
                fi
                SSRfilesize=`stat -c %s {inputSSR}`
                if [ $SSRfilesize = 0 ]
                then
                 touch {outputSSF} {outputSSR}
                fi
               Rscript ./scripts/match_pairs.r {inputASF},{inputASR},{inputSSF},{inputSSR} {outputASF},{outputASR},{outputSSF},{outputSSR} 
                 if grep -q "removed all reads: {outputASF}" ".gwf/logs/match_{project_name}_{library_id}_{tag_id}.stderr"
                 then
                  touch {outputASF}
                  touch {outputASR}
                 fi
                 if grep -q "removed all reads: {outputSSF}" ".gwf/logs/match_{project_name}_{library_id}_{tag_id}.stderr"
                 then                 
                  touch {outputSSF}
                  touch {outputSSR}
                 fi
            """.format(folderAS=folderAS,folderSS=folderSS,inputASF=input_files[0],inputASR=input_files[1],inputSSF=input_files[2],inputSSR=input_files[3],outputASF=output_files[0],outputASR=output_files[1],outputSSF=output_files[2],outputSSR=output_files[3],project_name=project_name,library_id=library_id,tag_id=tag_id)         

#Removing likely erroneous sequences 
for library_root in libraries:
    library_id = os.path.basename(library_root)
    
    input_files = []
    with open(os.path.join(library_root, "tags.txt")) as tags_file:
        for line in tags_file:
            tag_id, fseq, rseq = line.split()

            input_files.append("tmp/{}/DADA2_AS/filtered/matched/{}_F_matched.fastq.gz".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_AS/filtered/matched/{}_R_matched.fastq.gz".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_SS/filtered/matched/{}_F_matched.fastq.gz".format(library_id, tag_id))
            input_files.append("tmp/{}/DADA2_SS/filtered/matched/{}_R_matched.fastq.gz".format(library_id, tag_id))
            
    output_files = []
    
    output_files.append("tmp/{}/seqtab_AS_RDS".format(library_id))
    output_files.append("tmp/{}/seqtab.nochim_AS_RDS".format(library_id))
    output_files.append("tmp/{}/seqtab_SS_RDS".format(library_id))
    output_files.append("tmp/{}/seqtab.nochim_SS_RDS".format(library_id))
                                
    gwf.target(
      name="remove_errors_{}_{}".format(project_name,library_id),
      inputs=input_files,
      outputs=output_files,
      cores=4,
      memory="12g",
      walltime="4:00:00",
    ) << """
      Rscript ./scripts/remove_errors.r tmp/{library_id}
      """.format(library_id=library_id)         
            
#Summing sense and antisense sequence tables
for library_root in libraries:
    library_id = os.path.basename(library_root)
    input_files = []
    
    input_files.append("tmp/{}/seqtab_AS_RDS".format(library_id))
    input_files.append("tmp/{}/seqtab.nochim_AS_RDS".format(library_id))
    input_files.append("tmp/{}/seqtab_SS_RDS".format(library_id))
    input_files.append("tmp/{}/seqtab.nochim_SS_RDS".format(library_id))
    
    output_files = []
    
    output_files.append("tmp/{}/seqtab_RDS".format(library_id))
    output_files.append("tmp/{}/seqtab.nochim_RDS".format(library_id))
    
    gwf.target(
      name="sum_AS_SS_{}_{}".format(project_name,library_id),
      inputs=input_files,
      outputs=output_files,
      cores=1,
      memory="2g",
      walltime="4:00:00",
    ) << """
      Rscript ./scripts/sum_AS_SS.r tmp/{library_id}
      """.format(library_id=library_id)    

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
   cores=16,
   memory="32g",
   walltime="1:00:00",
 ) << """
   Rscript ./scripts/sum_libraries.r tmp/ results/
   """.format()                
