#!/usr/bin/env bash
set -euo pipefail

ANNOTATION_DIR=$(pwd)
EDTA_OUTPUT_DIR="edta_hapA"
REPEAT_OUTPUT_DIR="repeatmasker_abinitio"
GENOME_FILE="Vper.hapA.v1.unmasked.fasta"

# Caminhos de entrada
GENOME_PATH="$ANNOTATION_DIR/$EDTA_OUTPUT_DIR/$GENOME_FILE"
EDTA_GFF_PATH="$ANNOTATION_DIR/$EDTA_OUTPUT_DIR/${GENOME_FILE}.mod.EDTA.TEanno.gff3"

# Basename
GENOME_BASENAME="$(basename "$GENOME_FILE" .fasta)"
DB_NAME="${GENOME_BASENAME}_db"

THREADS=24

# Checagens
if ! command -v bedtools &> /dev/null; then
  echo "❌ Erro: BEDTools não encontrado. Ative o ambiente correto."
  exit 1
fi

if [ ! -s "$GENOME_PATH" ]; then
  echo "❌ Erro: FASTA não encontrado em $GENOME_PATH"
  exit 1
fi

mkdir -p "$REPEAT_OUTPUT_DIR"
cd "$REPEAT_OUTPUT_DIR"

echo "[1/3] Anotação de TEs com RepeatMasker e RepeatModeler..."

echo "-> BuildDatabase..."
BuildDatabase -name "$DB_NAME" "$GENOME_PATH"

echo "-> RepeatModeler..."
RepeatModeler -database "$DB_NAME" -threads "$THREADS" -LTRStruct || { echo "❌ RepeatModeler falhou."; exit 1; }

LIB="${DB_NAME}-families.fa"
if [ ! -s "$LIB" ]; then
  echo "❌ Erro: biblioteca do RepeatModeler não encontrada: $LIB"
  exit 1
fi

echo "-> RepeatMasker..."
# -pa para paralelismo; -dir para garantir saída aqui
RepeatMasker -pa "$THREADS" -lib "$LIB" -gff -dir "$REPEAT_OUTPUT_DIR" "$GENOME_PATH" || { echo "❌ RepeatMasker falhou."; exit 1; }

cd "$ANNOTATION_DIR"

# Caminho correto do GFF do RepeatMasker
REPEAT_GFF_PATH="$ANNOTATION_DIR/$REPEAT_OUTPUT_DIR/$REPEAT_OUTPUT_DIR/${GENOME_FILE}.out.gff"

echo "[2/3] Convertendo GFFs para BED..."

# RepeatMasker GFF -> BED (0-based), feature repeat_region
# Higieniza atributos para evitar espaços/; no campo 'name'
awk 'BEGIN{FS=OFS="\t"} 
     $0 !~ /^#/ && $3=="repeat_region" {
       name=$9; gsub(/ /,"_",name); gsub(/;/,"|",name);
       print $1, $4-1, $5, name, "1", $7
     }' "$REPEAT_GFF_PATH" > "$REPEAT_OUTPUT_DIR/repeatmasker_abinitio.bed" || { echo "❌ Falha ao converter GFF do RepeatMasker."; exit 1; }

# EDTA GFF -> BED
awk 'BEGIN{FS=OFS="\t"} 
     $0 !~ /^#/ && $3 ~ /transposable_element/ {
       name=$9; gsub(/ /,"_",name); gsub(/;/,"|",name);
       print $1, $4-1, $5, name, "1", $7
     }' "$EDTA_GFF_PATH" > "$EDTA_OUTPUT_DIR/edta_raw.bed" || { echo "❌ Falha ao converter GFF do EDTA."; exit 1; }

echo "[3/3] Mesclando anotações (não redundante) com BEDTools..."

# Ordena e funde (distância de 100 bp para colar hits próximos)
LC_ALL=C sort -k1,1 -k2,2n "$REPEAT_OUTPUT_DIR/repeatmasker_abinitio.bed" "$EDTA_OUTPUT_DIR/edta_raw.bed" \
| bedtools merge -i - -c 4,6 -o distinct,distinct -d 100 \
> "combined_te_annotations.bed" || { echo "❌ Falha na fusão com BEDTools."; exit 1; }

echo "✅ Concluído! Saída: $ANNOTATION_DIR/combined_te_annotations.bed"
