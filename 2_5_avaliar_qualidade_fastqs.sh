#!/bin/bash

# Diretório com os arquivos fastq combinados
INPUT_DIR=.
FASTQC_OUT=result_fastqc
MULTIQC_OUT=result_multiqc

# Cria as pastas de saída se não existirem
mkdir -p "$FASTQC_OUT"
mkdir -p "$MULTIQC_OUT"

# Lista todos os arquivos fastq.gz combinados
for fq in ${INPUT_DIR}/*_R[12].fastq.gz; do
    echo "Analisando qualidade de: $fq"
    fastqc "$fq" -o "$FASTQC_OUT"
done

# Rodar MultiQC sobre os resultados do FastQC
multiqc "$FASTQC_OUT" -o "$MULTIQC_OUT"

echo "✅ Análise de qualidade concluída."
echo "📁 Resultados do FastQC em: $FASTQC_OUT"
echo "📊 Sumário MultiQC em: $MULTIQC_OUT"
