#!/bin/bash

# Este script executa a anotação de TEs com RepeatModeler e RepeatMasker
# e combina os resultados com a anotação já gerada pelo EDTA.

# Diretórios e arquivos
# O diretório de trabalho principal onde os resultados finais serão salvos
ANNOTATION_DIR=$(pwd)
EDTA_OUTPUT_DIR="edta_hapA"
REPEAT_OUTPUT_DIR="repeatmasker_abinitio"
GENOME_FILE="Vper.hapA.v1.unmasked.fasta"

# Diretório de trabalho temporário em um disco rápido (ajuste se necessário)
# O `USER` é uma variável de ambiente que geralmente contém seu nome de usuário.
TMP_DIR="/scratch/$USER/te_annotation_$$" # O $$ garante um nome único

# Define o caminho completo para os arquivos
GENOME_PATH="$ANNOTATION_DIR/$EDTA_OUTPUT_DIR/$GENOME_FILE"
EDTA_GFF_PATH="$ANNOTATION_DIR/$EDTA_OUTPUT_DIR/${GENOME_FILE}.mod.EDTA.TEanno.gff3"

GENOME_BASENAME=$(basename $GENOME_FILE .fasta)
THREADS=24

# Verifique se o BEDTools está instalado
if ! command -v bedtools &> /dev/null; then
    echo "❌ Erro: BEDTools não encontrado. Por favor, instale-o ou ative o ambiente Conda."
    exit 1
fi

# Cria o diretório temporário e copia o genoma para lá
echo "Criando diretório temporário em $TMP_DIR e copiando o genoma..."
mkdir -p "$TMP_DIR" || { echo "❌ Erro: Não foi possível criar o diretório temporário. Verifique as permissões."; exit 1; }
cp "$GENOME_PATH" "$TMP_DIR/" || { echo "❌ Erro: Não foi possível copiar o genoma para o diretório temporário."; exit 1; }

# Muda para o diretório temporário para otimizar o I/O
echo "Mudando para o diretório temporário para a execução..."
cd "$TMP_DIR" || { echo "❌ Erro: Não foi possível mudar para o diretório temporário."; exit 1; }

# =========================================================
# Passo 1: Anotação de Homologia e ab initio com RepeatMasker/RepeatModeler
# =========================================================
echo "[1/3] Anotação de TEs com RepeatMasker e RepeatModeler..."

# Rodar RepeatModeler para criar uma biblioteca de repetições
echo "-> Rodando RepeatModeler..."
BuildDatabase -name "${GENOME_BASENAME}_db" "$GENOME_FILE" || { echo "❌ Erro: BuildDatabase falhou."; exit 1; }
RepeatModeler -database "${GENOME_BASENAME}_db" -threads "$THREADS" -LTRStruct || { echo "❌ Erro: RepeatModeler falhou."; exit 1; }

# Rodar RepeatMasker com a biblioteca customizada
echo "-> Rodando RepeatMasker com a biblioteca gerada..."
RepeatMasker -threads "$THREADS" -lib "${GENOME_BASENAME}_db-families.fa" -gff "$GENOME_FILE" || { echo "❌ Erro: RepeatMasker falhou."; exit 1; }

# Voltando para o diretório original
echo "Movendo resultados e limpando o diretório temporário..."
cd "$ANNOTATION_DIR" || { echo "❌ Erro: Não foi possível voltar para o diretório de trabalho original."; exit 1; }

# Cria o diretório de saída e move os resultados para lá
mkdir -p "$REPEAT_OUTPUT_DIR"
mv "$TMP_DIR/${GENOME_FILE}.gff" "$REPEAT_OUTPUT_DIR/"
mv "$TMP_DIR/${GENOME_BASENAME}_db-families.fa" "$REPEAT_OUTPUT_DIR/"
mv "$TMP_DIR/${GENOME_FILE}.tbl" "$REPEAT_OUTPUT_DIR/"

# Removendo o diretório temporário
rm -rf "$TMP_DIR"

# =========================================================
# Passo 2: Conversão e Fusão das Anotações
# =========================================================
echo "[2/3] Mesclando anotações de forma não redundante com BEDTools..."

# Define o caminho para o GFF gerado na etapa anterior
REPEAT_GFF_PATH="$ANNOTATION_DIR/$REPEAT_OUTPUT_DIR/${GENOME_FILE}.gff"

# Converte o GFF do RepeatMasker para BED
echo "-> Convertendo anotação do RepeatMasker para BED..."
grep -v "^#" "$REPEAT_GFF_PATH" | awk 'BEGIN{FS="\t"; OFS="\t"}{if($3=="repeat") print $1, $4-1, $5, $9, "1", $7}' > "$REPEAT_OUTPUT_DIR/repeatmasker_abinitio.bed" || { echo "❌ Erro: Falha ao converter GFF do RepeatMasker."; exit 1; }

# Converte a anotação do EDTA para BED
echo "-> Convertendo anotação do EDTA para BED..."
grep -v "^#" "$EDTA_GFF_PATH" | awk 'BEGIN{FS="\t"; OFS="\t"}{if($3 ~ /transposable_element/) print $1, $4-1, $5, $9, "1", $7}' > "$EDTA_OUTPUT_DIR/edta_raw.bed" || { echo "❌ Erro: Falha ao converter GFF do EDTA."; exit 1; }

# Combina as duas anotações e resolve as sobreposições com BEDTools
echo "-> Realizando fusão com bedtools..."
cat "$REPEAT_OUTPUT_DIR/repeatmasker_abinitio.bed" "$EDTA_OUTPUT_DIR/edta_raw.bed" | sort -k1,1 -k2,2n | bedtools merge -i - -c 4,6 -o distinct,distinct -d 100 > "combined_te_annotations.bed" || { echo "❌ Erro: Falha na fusão com BEDTools."; exit 1; }

echo "✅ Etapa de fusão concluída! O arquivo final é 'combined_te_annotations.bed'."

# =========================================================
# Passo 3: Identificar o genoma "soft-masked"
# =========================================================
echo "[3/3] Localizando o genoma 'soft-masked' para o BRAKER3..."

SOFT_MASKED_GENOME="$EDTA_OUTPUT_DIR/${GENOME_BASENAME}.mod.MAKER.masked"

echo "✅ O arquivo FASTA do genoma que deve ser usado como input para o BRAKER3 é:"
echo "$SOFT_MASKED_GENOME"
