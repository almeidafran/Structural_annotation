#!/bin/bash

# Caminho onde estão os arquivos de RNA-seq
INPUT_DIR=~/rna_seq_peripherica/raw_data

# Lista de prefixos das amostras
samples=(VP1C VP1F VP1R VP2C VP2F VP2R VP3C VP3F VP3R)

# Loop para combinar os arquivos
for sample in "${samples[@]}"; do
    echo "Combinando arquivos da amostra: $sample"

    # Combinando R1
    cat ${INPUT_DIR}/${sample}_S*_L001_R1_001.fastq.gz ${INPUT_DIR}/${sample}_S*_L002_R1_001.fastq.gz > ${sample}_R1.fastq.gz

    # Combinando R2
    cat ${INPUT_DIR}/${sample}_S*_L001_R2_001.fastq.gz ${INPUT_DIR}/${sample}_S*_L002_R2_001.fastq.gz > ${sample}_R2.fastq.gz

    echo "Gerado: ${sample}_R1.fastq.gz e ${sample}_R2.fastq.gz"
    echo "------------------------------------------------------"
done
