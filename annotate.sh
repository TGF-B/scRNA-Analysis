#!/bin/bash
# filepath: /media/donaldtangai4s/Yu_Omics2/fastq/annotate.sh
cellranger="/home/donaldtangai4s/cellranger_temp/cellranger-9.0.0/bin/cellranger"

# List of samples
samples=("Ctrl-3")

# Annotate each sample
for sample in "${samples[@]}"; do
    echo "Annotating cellranger output for sample: $sample"
    $cellranger annotate \
        --id="$sample" \
        --matrix="/home/donaldtangai4s/new/$sample/outs/filtered_feature_bc_matrix.h5" \
        --cell-annotation-model=mouse_pca_v1_beta \
        --tenx-cloud-token-path="/media/donaldtangai4s/Yu_Omics2/fastq/token.txt" \
        --output-dir="/home/donaldtangai4s/new/$sample/annotated" \
        --localcores=10 \
        --localmem=5

echo "Annotation completed for all samples."
done