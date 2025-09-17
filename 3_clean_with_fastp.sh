#!/bin/bash

# Diretórios
INPUT_DIR=~/rna_seq_peripherica/merged_data
OUTPUT_DIR=~/rna_seq_peripherica/cleaned_fastp
mkdir -p "$OUTPUT_DIR"

# Loop para todos os arquivos R1
for R1 in "$INPUT_DIR"/*_R1.fastq.gz; do
  # Derivar o nome base da amostra
  SAMPLE=$(basename "$R1" _R1.fastq.gz)
  R2="$INPUT_DIR/${SAMPLE}_R2.fastq.gz"

  # Arquivos de saída
  OUT_R1="$OUTPUT_DIR/${SAMPLE}_R1_clean.fq.gz"
  OUT_R2="$OUTPUT_DIR/${SAMPLE}_R2_clean.fq.gz"
  HTML_REPORT="$OUTPUT_DIR/${SAMPLE}_fastp.html"
  JSON_REPORT="$OUTPUT_DIR/${SAMPLE}_fastp.json"

  echo "🔄 Processando $SAMPLE ..."

  # Rodar fastp
  fastp \
    -i "$R1" -I "$R2" \
    -o "$OUT_R1" -O "$OUT_R2" \
    --detect_adapter_for_pe \
    --thread 8 \
    --n_base_limit 5 \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --html "$HTML_REPORT" \
    --json "$JSON_REPORT"

  echo "✅ Finalizado: $SAMPLE"
done

echo "🎉 Todos os arquivos foram processados com fastp."
