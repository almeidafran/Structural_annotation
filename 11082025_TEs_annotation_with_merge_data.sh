#!/bin/bash

# Diretório de trabalho e arquivos
GENOME="22042025_renomeado_hap1_only_9chr_peripherica.fa"
GENOME_BASENAME=$(basename $GENOME .fa)
THREADS=24

# Verifique se o BEDTools está instalado e no seu PATH
if ! command -v bedtools &> /dev/null; then
    echo "❌ Erro: BEDTools não encontrado. Por favor, instale-o ou ative o ambiente Conda."
    exit 1
fi

# =========================================================
# Passo 1: Anotação de Homologia e ab initio com RepeatMasker/RepeatModeler
# =========================================================
echo "[1/4] Anotação de TEs com RepeatMasker e RepeatModeler..."

# Rodar RepeatModeler para criar uma biblioteca de repetições
echo "-> Rodando RepeatModeler..."
BuildDatabase -name ${GENOME_BASENAME}_db $GENOME
RepeatModeler -database ${GENOME_BASENAME}_db -pa $THREADS -LTRStruct

# Rodar RepeatMasker com a biblioteca customizada
echo "-> Rodando RepeatMasker com a biblioteca gerada..."
RepeatMasker -pa $THREADS -lib ${GENOME_BASENAME}_db-families.fa -gff -dir repeatmasker_abinitio $GENOME

# =========================================================
# Passo 2: Anotação com EDTA
# =========================================================
echo "[2/4] Anotação de TEs com EDTA..."
EDTA.pl \
    --genome $GENOME \
    --species others \
    --step all \
    --threads $THREADS \
    --overwrite 0 \
    --sensitive 1 \
    --anno 1 \
    --evaluate 0 \
    --cds NA \
    --rna NA

# =========================================================
# Passo 3: Conversão e Fusão das Anotações
# =========================================================
echo "[3/4] Mesclando anotações de forma não redundante com BEDTools..."

# Converte o GFF do RepeatMasker para BED
# A saída do RepeatMasker em GFF precisa ser tratada
echo "-> Convertendo anotação do RepeatMasker para BED..."
gff_repeatmasker="repeatmasker_abinitio/${GENOME_BASENAME}.fa.gff"
grep -v "^#" "$gff_repeatmasker" | awk 'BEGIN{FS="\t"; OFS="\t"}{if($3=="repeat") print $1, $4-1, $5, $9, "1", $7}' > repeatmasker_abinitio.bed

# O EDTA gera um arquivo GFF completo na pasta de saída
echo "-> Convertendo anotação do EDTA para BED..."
gff_edta="EDTA_raw/${GENOME_BASENAME}.mod.EDTA.TE.gff3"
grep -v "^#" "$gff_edta" | awk 'BEGIN{FS="\t"; OFS="\t"}{if($3 ~ /transposable_element/) print $1, $4-1, $5, $9, "1", $7}' > edta_raw.bed

# Combina as duas anotações e resolve as sobreposições.
# Prioridade para o EDTA, pois é geralmente mais completo.
# A ferramenta 'bedtools' merge junta sobreposições e unifica as anotações
echo "-> Realizando fusão com bedtools..."
cat repeatmasker_abinitio.bed edta_raw.bed | sort -k1,1 -k2,2n | bedtools merge -i - -c 4,6 -o distinct,distinct -d 100 > combined_te_annotations.bed

echo "✅ Etapa de fusão concluída! O arquivo final é 'combined_te_annotations.bed'."

# =========================================================
# Passo 4: Limpeza (Opcional)
# =========================================================
echo "[4/4] Limpando arquivos intermediários..."
# rm -f repeatmasker_abinitio.bed edta_raw.bed
echo "✅ Finalizado! O arquivo final de anotação de TEs está em combined_te_annotations.bed."
