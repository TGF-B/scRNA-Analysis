#!/bin/bash
# chmod +x cellranger.sh   /before running below
# transcriptome 
transcriptome_path="media/i/Yu_Omics2/refdata-gex-mm10-2020-A"

# fastq 
fastq_path="media/i/Yu_Omics2"

# fastq
fastq_files=(
    "CP-1_S1_L001_R1_001.fastq"
    "CP-1_S1_L001_R2_001.fastq"
    "CP-2_S1_L001_R1_001.fastq"
    "CP-2_S1_L001_R2_001.fastq"
    "CP-3_S1_L001_R1_001.fastq"
    "CP-3_S1_L001_R2_001.fastq"
    "Ctrl-1_S1_L001_R1_001.fastq"
    "Ctrl-1_S1_L001_R2_001.fastq"
    "Ctrl-2_S1_L001_R1_001.fastq"
    "Ctrl-2_S1_L001_R2_001.fastq"
    "Ctrl-3_S1_L001_R1_001.fastq"
    "Ctrl-3_S1_L001_R2_001.fastq"
)

# tranverse all the fastqs
for fastq_file in "${fastq_files[@]}"; do
    # name id and sample with former 5 letters of fastq`s filename
    id_sample=${fastq_file:0:5}

    # executing Cell Ranger
    cellranger count --id="$id_sample" --transcriptome="$transcriptome_path" --fastqs="$fastq_path" --sample="$id_sample"

done

#single command
# cellranger count --id=Ctrl3a --transcriptome=./refdata-gex-mm10-2020-A --fastqs=./ --sample=Ctrl-2b
