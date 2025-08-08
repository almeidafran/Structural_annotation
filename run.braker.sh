#!/bin/bash

# Caminhos
WORKDIR=~/rna_seq_peripherica

# Caminho para o genoma
GENOME=$(dirname "$WORKDIR")/2306_assembly_Todd/Vper.hapA.v1.munged/Vper.hapA.v1.softmasked.fasta

# Diretório com os arquivos BAM gerados pelo HISAT2
ALIGN_DIR=$WORKDIR/results_hisat2_genoma_todd_hapA

# Diretório de saída para o BRAKER3
BRAKER_DIR=$WORKDIR/results_braker3_hapA
mkdir -p "$BRAKER_DIR"

# Caminho para evidências de proteínas (OPCIONAL)
# Se você tiver um arquivo fasta de proteínas de uma espécie relacionada,
# substitua o caminho abaixo. Caso contrário, deixe como está (caminho vazio).
PROT_SEQ=""

# Parâmetros de execução
THREADS=8  # Número de threads a serem usadas. Ajuste conforme a sua necessidade.
SPECIES_NAME="Vper_hapA" # Nome da espécie para o BRAKER.

# Coleta todos os arquivos BAM do diretório de alinhamento
BAMS=$(find "$ALIGN_DIR" -name "*.sorted.bam" | tr '\n' ',' | sed 's/,$//')

# Verifique se há arquivos BAM
if [ -z "$BAMS" ]; then
    echo "❌ Erro: Nenhum arquivo .sorted.bam encontrado em $ALIGN_DIR"
    exit 1
fi

echo "### Etapa 1: Predição de genes com BRAKER3..."
echo "✅ Arquivos BAM encontrados: $BAMS"
echo "✅ Diretório de saída do BRAKER: $BRAKER_DIR"

# Executa o comando BRAKER3
# --genome: caminho para o genoma
# --hints: evidências de RNA-seq (arquivos BAM)
# --softmasking: indica que o genoma já está soft-masked
# --species: nome da espécie
# --threads: número de threads para acelerar a execução
# --workingdir: diretório onde os resultados serão salvos

if [ -z "$PROT_SEQ" ]; then
    echo "⚠️ Aviso: Nenhuma evidência de proteína foi fornecida. Rodando BRAKER3 apenas com dados de RNA-seq."
    braker.pl \
        --genome="$GENOME" \
        --hints="$BAMS" \
        --softmasking \
        --species="$SPECIES_NAME" \
        --threads "$THREADS" \
        --workingdir="$BRAKER_DIR"
else
    echo "✅ Usando evidências de proteína do arquivo: $PROT_SEQ"
    braker.pl \
        --genome="$GENOME" \
        --hints="$BAMS" \
        --prot_seq="$PROT_SEQ" \
        --softmasking \
        --species="$SPECIES_NAME" \
        --threads "$THREADS" \
        --workingdir="$BRAKER_DIR"
fi

echo "### Pipeline do BRAKER3 finalizada com sucesso!"
