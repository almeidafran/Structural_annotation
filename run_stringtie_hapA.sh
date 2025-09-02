#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Uso:
  run_stringtie_hapA.sh \
    -b <dir_bams> \
    -g <genoma.fasta> \
    -o <outdir> \
    -t <threads> \
    -s <FR|RF|UNSTRANDED> \
    [--min-cov N] [--min-frac X] [--min-junc N] [--min-len N]

Obrigatórios:
  -b  Diretório contendo arquivos *.sorted.bam (e seus .bai)
  -g  Caminho para o FASTA do genoma
  -o  Diretório de saída
  -t  Número de threads
  -s  Tipo de biblioteca (FR, RF ou UNSTRANDED)

O que faz:
  (1) Roda StringTie em cada *.sorted.bam do diretório
  (2) Faz o merge dos GTFs (stringtie --merge)
  (3) Extrai cDNAs do genoma com gffread

Parâmetros finos (opcionais; defaults entre colchetes):
  --min-cov N     (-c) cobertura mínima por 'bundle'          [2]
  --min-frac X    (-f) fração mínima de isoforma (TPM)        [0.1]
  --min-junc N    (-j) cobertura mínima em junções            [2]
  --min-len N     (-m) tamanho mínimo de transcrito (nt)      [200]

Exemplo:
  ./run_stringtie_hapA.sh \
    -b ~/rna_seq_peripherica/results_hisat2_log_genoma_todd_hapA \
    -g /dados/genomas/Vper.hapA.v1.fasta \
    -o ~/rna_seq_peripherica/stringtie_hapA_out \
    -t 16 \
    -s FR \
    --min-cov 2 --min-frac 0.1 --min-junc 2 --min-len 200
USAGE
}

# defaults dos parâmetros finos
MIN_COVERAGE=2
MIN_FRAC=0.1
MIN_JUNC_COV=2
MIN_LENGTH=200

# lê flags curtas com getopts e flags longas manualmente
BAM_DIR=""
GENOME_FA=""
OUTDIR=""
THREADS=""
LIBTYPE=""

# pré-processa flags longas → curtas auxiliares
LONGOPTS=(
  --min-cov
  --min-frac
  --min-junc
  --min-len
  --help
)
# varre argumentos para capturar longos
PARSED_ARGS=()
while (( "$#" )); do
  case "$1" in
    --min-cov)   shift; MIN_COVERAGE="${1:-}";;
    --min-frac)  shift; MIN_FRAC="${1:-}";;
    --min-junc)  shift; MIN_JUNC_COV="${1:-}";;
    --min-len)   shift; MIN_LENGTH="${1:-}";;
    --help) usage; exit 0;;
    -*)
      PARSED_ARGS+=("$1" "${2:-}")
      shift
      ;;
    *)
      PARSED_ARGS+=("$1")
      ;;
  esac
  shift || true
done

# repassa os restantes (curtos) para getopts
set -- "${PARSED_ARGS[@]}"

while getopts ":b:g:o:t:s:h" opt; do
  case "${opt}" in
    b) BAM_DIR="${OPTARG}";;
    g) GENOME_FA="${OPTARG}";;
    o) OUTDIR="${OPTARG}";;
    t) THREADS="${OPTARG}";;
    s) LIBTYPE="${OPTARG}";;
    h) usage; exit 0;;
    \?) echo "❌ Opção inválida: -${OPTARG}"; usage; exit 2;;
    :)  echo "❌ Faltando argumento para -${OPTARG}"; usage; exit 2;;
  esac
done

# valida obrigatórios
[[ -z "${BAM_DIR}"   ]] && echo "❌ Falta -b <dir_bams>" && usage && exit 2
[[ -z "${GENOME_FA}" ]] && echo "❌ Falta -g <genoma.fasta>" && usage && exit 2
[[ -z "${OUTDIR}"    ]] && echo "❌ Falta -o <outdir>" && usage && exit 2
[[ -z "${THREADS}"   ]] && echo "❌ Falta -t <threads>" && usage && exit 2
[[ -z "${LIBTYPE}"   ]] && echo "❌ Falta -s <FR|RF|UNSTRANDED>" && usage && exit 2

# checagens de binários
command -v stringtie >/dev/null 2>&1 || { echo "❌ stringtie não encontrado no PATH"; exit 1; }
command -v gffread   >/dev/null 2>&1 || { echo "❌ gffread não encontrado no PATH"; exit 1; }
command -v samtools  >/dev/null 2>&1 || { echo "❌ samtools não encontrado no PATH"; exit 1; }

# checagens de caminhos
[[ -d "${BAM_DIR}"    ]] || { echo "❌ BAM_DIR não existe: ${BAM_DIR}"; exit 1; }
[[ -f "${GENOME_FA}"  ]] || { echo "❌ GENOME_FA não existe: ${GENOME_FA}"; exit 1; }

mkdir -p "${OUTDIR}/per_sample" "${OUTDIR}/logs"

# mapeia strand flag
declare -a STRAND_FLAG=()
case "${LIBTYPE^^}" in
  FR) STRAND_FLAG=(--fr) ;;
  RF) STRAND_FLAG=(--rf) ;;
  UNSTRANDED) STRAND_FLAG=() ;;
  *) echo "⚠️  LIBTYPE='${LIBTYPE}' inválido. Use FR, RF ou UNSTRANDED. Continuando sem strand."; STRAND_FLAG=() ;;
esac

echo "==> Config:"
echo "    BAM_DIR       : ${BAM_DIR}"
echo "    GENOME_FA     : ${GENOME_FA}"
echo "    OUTDIR        : ${OUTDIR}"
echo "    THREADS       : ${THREADS}"
echo "    LIBTYPE       : ${LIBTYPE}"
echo "    MIN_COVERAGE  : ${MIN_COVERAGE}"
echo "    MIN_FRAC      : ${MIN_FRAC}"
echo "    MIN_JUNC_COV  : ${MIN_JUNC_COV}"
echo "    MIN_LENGTH    : ${MIN_LENGTH}"

MERGELIST="${OUTDIR}/mergelist.txt"
: > "${MERGELIST}"

shopt -s nullglob
BAMS=( "${BAM_DIR}"/*.sorted.bam )
shopt -u nullglob
(( ${#BAMS[@]} )) || { echo "❌ Nenhum *.sorted.bam em ${BAM_DIR}"; exit 1; }

# 1) StringTie por amostra
for BAM in "${BAMS[@]}"; do
  SAMPLE="$(basename "${BAM}" .sorted.bam)"
  GTF_OUT="${OUTDIR}/per_sample/${SAMPLE}.gtf"
  LOG_OUT="${OUTDIR}/logs/${SAMPLE}.log"

  # garante índice
  if [[ ! -f "${BAM}.bai" ]]; then
    echo "   - Indexando ${BAM}..."
    samtools index -@ "${THREADS}" "${BAM}"
  fi

  echo ">> StringTie: ${SAMPLE}"
  {
    echo "[INFO] $(date) - StringTie ${SAMPLE}"
    echo "[CMD ] stringtie ${BAM} -p ${THREADS} -o ${GTF_OUT} ${STRAND_FLAG[*]} -c ${MIN_COVERAGE} -f ${MIN_FRAC} -j ${MIN_JUNC_COV} -m ${MIN_LENGTH}"
  } > "${LOG_OUT}"

  stringtie "${BAM}" \
    -p "${THREADS}" \
    -o "${GTF_OUT}" \
    "${STRAND_FLAG[@]}" \
    -c "${MIN_COVERAGE}" \
    -f "${MIN_FRAC}" \
    -j "${MIN_JUNC_COV}" \
    -m "${MIN_LENGTH}" \
    >> "${LOG_OUT}" 2>&1

  echo "${GTF_OUT}" >> "${MERGELIST}"
done

# 2) Merge GTFs
MERGED_GTF="${OUTDIR}/stringtie_merged.gtf"
echo ">> Merge dos GTFs → ${MERGED_GTF}"
stringtie --merge -p "${THREADS}" -o "${MERGED_GTF}" "${MERGELIST}" \
  > "${OUTDIR}/logs/merge.log" 2>&1

[[ -f "${MERGED_GTF}" ]] || { echo "❌ Falhou em gerar ${MERGED_GTF}"; exit 1; }

# 3) Extrai cDNAs
TRANSCRIPTS_FA="${OUTDIR}/stringtie_transcripts.fa"
echo ">> gffread → ${TRANSCRIPTS_FA}"
gffread "${MERGED_GTF}" -g "${GENOME_FA}" -w "${TRANSCRIPTS_FA}"

echo "✅ Concluído!"
echo "   - GTFs por amostra: ${OUTDIR}/per_sample/*.gtf"
echo "   - Merge            : ${MERGED_GTF}"
echo "   - Transcritos (FA) : ${TRANSCRIPTS_FA}"
