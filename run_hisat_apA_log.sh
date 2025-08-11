#!/bin/bash

# Caminhos
WORKDIR=$(pwd)
TRIMMED_DIR=$WORKDIR/cleaned_40_fastp

# Caminho completo para o genoma
GENOME=$(dirname "$WORKDIR")/2306_assembly_Todd/Vper.hapA.v1.munged/Vper.hapA.v1.softmasked.fasta

# Caminho para o índice
GENOME_INDEX=$(dirname "$WORKDIR")/2306_assembly_Todd/Vper.hapA.v1.munged/Vper.hapA.v1.munged

# Diretórios de saída
ALIGN_DIR=$WORKDIR/results_hisat2_log_genoma_todd_hapA
LOG_DIR=$ALIGN_DIR/logs

mkdir -p "$ALIGN_DIR" "$LOG_DIR"

# 1. Criar índice HISAT2 (se ainda não existir)
echo "### Etapa 1: Criar índice do HISAT2..."
if [ ! -f "${GENOME_INDEX}.1.ht2" ]; then
    hisat2-build "$GENOME" "$GENOME_INDEX"
fi

# 2. Alinhamento com HISAT2 + sort + index
echo "### Etapa 2: Alinhamento com HISAT2..."
cd "$TRIMMED_DIR"

for R1 in *_R1_clean.fq.gz; do
    BASE=$(basename "$R1" _R1_clean.fq.gz)
    R2="${BASE}_R2_clean.fq.gz"
    OUT_BAM="$ALIGN_DIR/${BASE}.sorted.bam"
    LOG_FILE="$LOG_DIR/${BASE}.log"

    echo "Alinhando $BASE..."
    
    hisat2 -x "$GENOME_INDEX" -1 "$R1" -2 "$R2" \
        -p 4 --dta 2> "$LOG_FILE" | \
        samtools sort -@ 4 -o "$OUT_BAM"

    samtools index "$OUT_BAM"
done

echo "### Pipeline finalizada com sucesso!"
