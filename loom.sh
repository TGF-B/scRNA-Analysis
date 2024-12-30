#!/bin/bash

# List of samples
samples=("CP-2" "CP-3" "Ctrl-1" "Ctrl-2" "Ctrl-3")
# Output directory for scVelo results
scvelo_output_dir="/home/donaldtangai4s/new"
mkdir -p "$scvelo_output_dir"

# Path to the annotation GTF file
gtf_file="/media/donaldtangai4s/Yu_Omics2/fastq/refdata-gex-mm10-2020-A/genes/genes.gtf"
velocyto="/home/donaldtangai4s/miniconda3/envs/py39/bin/velocyto"
python="/home/donaldtangai4s/miniconda3/envs/py39/bin/python3" #python 3.9.20
# Run Velocyto for each sample
for sample in "${samples[@]}"; do
    sample_output_dir="/home/donaldtangai4s/new/$sample/outs"
    sample_scvelo_output_dir="$scvelo_output_dir/$sample"
    mkdir -p "$sample_scvelo_output_dir"
    
    echo "Running Velocyto for sample: $sample"
    $velocyto run -b "$sample_output_dir/filtered_feature_bc_matrix/barcodes.tsv.gz" \
                  -o "$sample_output_dir" \
                    "$sample_output_dir/possorted_genome_bam.bam" \
                    "$gtf_file"
    
    loom_file="$sample_output_dir/${sample}.loom"
    
    echo "scVelo analysis completed for sample: $sample"
done



#Velocyto parsing the wrong python path