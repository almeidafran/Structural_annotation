#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Pipeline de QC em 3 fases (pré-PASA, pós-PASA, final EVM)

Uso:
  qc_all_phases.sh \
    --genome <genome.fa> \
    --bamdir <bam_dir> \
    --noG <st_noG.gtf> \
    --stB <st_braker.gtf> \
    --stBH <st_braker_helixer.gtf> \
    --pasa <pasa_refined.gff3> \
    --stChosen <stringtie_escolhido.gtf> \
    --evm <evm_final.gff3> \
    --lineage <busco_lineage> \
    -t <threads> \
    -o <outdir_raiz>

Cria:
  <outdir_raiz>/01_prePASA/
  <outdir_raiz>/02_posPASA/
  <outdir_raiz>/03_finalEVM/
USAGE
}

GENOME="" ; BAMDIR="" ; NO_GTF="" ; STB_GTF="" ; STBH_GTF="" ; PASA_GFF="" ; ST_CHOSEN="" ; EVM_GFF="" ; LINEAGE="" ; THREADS=8 ; ROOT=""
while (( "$#" )); do
  case "$1" in
    --genome) shift; GENOME="${1:-}";;
    --bamdir) shift; BAMDIR="${1:-}";;
    --noG) shift; NO_GTF="${1:-}";;
    --stB) shift; STB_GTF="${1:-}";;
    --stBH) shift; STBH_GTF="${1:-}";;
    --pasa) shift; PASA_GFF="${1:-}";;
    --stChosen) shift; ST_CHOSEN="${1:-}";;
    --evm) shift; EVM_GFF="${1:-}";;
    --lineage) shift; LINEAGE="${1:-}";;
    -t) shift; THREADS="${1:-}";;
    -o) shift; ROOT="${1:-}";;
    -h|--help) usage; exit 0;;
    *) echo "arg desconhecido: $1"; usage; exit 2;;
  esac; shift || true
done

[[ -f "$GENOME" && -d "$BAMDIR" && -f "$NO_GTF" && -f "$STB_GTF" && -f "$STBH_GTF" && -f "$PASA_GFF" && -f "$ST_CHOSEN" && -f "$EVM_GFF" && -n "$LINEAGE" && -n "$ROOT" ]] || { usage; exit 2; }

PRE="$ROOT/01_prePASA"
POS="$ROOT/02_posPASA"
FIN="$ROOT/03_finalEVM"
mkdir -p "$PRE" "$POS" "$FIN"

# 1) Pre-PASA
01_prePASA_compare_stringtie.sh -g "$GENOME" -b "$BAMDIR" --noG "$NO_GTF" --stB "$STB_GTF" --stBH "$STBH_GTF" -t "$THREADS" -o "$PRE" --lineage "$LINEAGE"

# 2) Pós-PASA
02_posPASA_qc.sh -g "$GENOME" -a "$PASA_GFF" -x "$ST_CHOSEN" -b "$BAMDIR" -t "$THREADS" -o "$POS"

# 3) Final EVM
03_finalEVM_qc.sh -g "$GENOME" -a "$EVM_GFF" -b "$BAMDIR" -t "$THREADS" -l "$LINEAGE" -o "$FIN"

echo "🎉 Tudo pronto! Pastas: $PRE | $POS | $FIN"
