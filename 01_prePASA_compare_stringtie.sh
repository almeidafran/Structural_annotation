#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
[01] Pré-PASA: comparar 3 GTFs do StringTie (noG, BRAKER, BRAKER+Helixer)

Uso:
  01_prePASA_compare_stringtie.sh \
    -g <genome.fa> \
    -b <bam_dir> \
    --noG <st_noG.gtf> \
    --stB <st_braker.gtf> \
    --stBH <st_braker_helixer.gtf> \
    -t <threads> \
    -o <outdir> \
    [--lineage embryophyta_odb10]  # opcional: BUSCO leve (transcriptome)

Saídas:
  <outdir>/
    01_stats/agat_*.txt                 (estatísticas estruturais)
    02_busco/ (se --lineage)           (BUSCO transcriptome)
    03_junctions/summary.txt           (% introns suportados por junções RNA-seq)
    04_exon_support/summary.txt        (% exons com cov > 0 por RNA-seq)
USAGE
}

GENOME="" ; BAM_DIR="" ; NO_GTF="" ; STB_GTF="" ; STBH_GTF="" ; THREADS=8 ; OUTDIR="" ; LINEAGE=""
PARSED=()
while (( "$#" )); do
  case "$1" in
    --noG)  shift; NO_GTF="${1:-}";;
    --stB)  shift; STB_GTF="${1:-}";;
    --stBH) shift; STBH_GTF="${1:-}";;
    --lineage) shift; LINEAGE="${1:-}";;
    -g) shift; GENOME="${1:-}";;
    -b) shift; BAM_DIR="${1:-}";;
    -t) shift; THREADS="${1:-}";;
    -o) shift; OUTDIR="${1:-}";;
    -h|--help) usage; exit 0;;
    *) PARSED+=("$1");;
  esac; shift || true
done
set -- "${PARSED[@]}"

[[ -f "$GENOME" && -d "$BAM_DIR" && -f "$NO_GTF" && -f "$STB_GTF" && -f "$STBH_GTF" && -n "$OUTDIR" ]] || { usage; exit 2; }

for bin in agat_sp_statistics.pl gffread bedtools regtools gffcompare; do
  command -v "$bin" >/dev/null 2>&1 || { echo "❌ falta $bin no PATH"; exit 1; }
done
[[ -n "$LINEAGE" ]] && command -v busco >/dev/null 2>&1 || true

mkdir -p "$OUTDIR"/{01_stats,02_busco,03_junctions/junctions_bed12,04_exon_support,logs}

# 1) AGAT stats
declare -A MAP=( ["noG"]="$NO_GTF" ["stB"]="$STB_GTF" ["stBH"]="$STBH_GTF" )
for tag in noG stB stBH; do
  GTF="${MAP[$tag]}"
  agat_sp_statistics.pl --gff "$GTF" -g "$GENOME" \
    -o "$OUTDIR/01_stats/agat_${tag}.txt" > "$OUTDIR/logs/agat_${tag}.log" 2>&1
done

# 2) BUSCO (opcional - transcriptome) a partir de transcritos extraídos
if [[ -n "$LINEAGE" ]]; then
  for tag in noG stB stBH; do
    GTF="${MAP[$tag]}"
    gffread "$GTF" -g "$GENOME" -w "$OUTDIR/02_busco/${tag}.transcripts.fa"
    busco -i "$OUTDIR/02_busco/${tag}.transcripts.fa" -l "$LINEAGE" -m transcriptome -c "$THREADS" \
      -o "busco_tx_${tag}" --out_path "$OUTDIR/02_busco" > "$OUTDIR/logs/busco_tx_${tag}.log" 2>&1 || true
  done
fi

# 3) Suporte de junções (RNA-seq) → % introns suportados
for BAM in "$BAM_DIR"/*.sorted.bam; do
  S=$(basename "$BAM" .sorted.bam)
  regtools junctions extract -a 8 -m 50 -M 500000 "$BAM" > "$OUTDIR/03_junctions/junctions_bed12/${S}.bed12"
done
cat "$OUTDIR"/03_junctions/junctions_bed12/*.bed12 > "$OUTDIR/03_junctions/all_junctions.bed12"

for tag in noG stB stBH; do
  GTF="${MAP[$tag]}"
  # introns a partir do GTF
  agat_convert_sp_gff2bed.pl --gff "$GTF" --out "$OUTDIR/03_junctions/${tag}.bed" > /dev/null 2>&1
  awk '$8=="intron"{OFS="\t"; print $1,$2,$3,$4,$5,$6}' "$OUTDIR/03_junctions/${tag}.bed" > "$OUTDIR/03_junctions/${tag}.introns.bed"
  bedtools intersect -u -a "$OUTDIR/03_junctions/${tag}.introns.bed" -b "$OUTDIR/03_junctions/all_junctions.bed12" \
    | wc -l > "$OUTDIR/03_junctions/introns_supported_${tag}.count"
  wc -l "$OUTDIR/03_junctions/${tag}.introns.bed" | awk '{print $1}' > "$OUTDIR/03_junctions/introns_total_${tag}.count"
done

python - <<'PY'
import os
base=os.environ["BASE"]
tags=["noG","stB","stBH"]
with open(os.path.join(base,"summary.txt"),"w") as out:
  for tag in tags:
    sup=int(open(os.path.join(base,f"introns_supported_{tag}.count")).read().strip())
    tot=int(open(os.path.join(base,f"introns_total_{tag}.count")).read().strip())
    pct=100.0*sup/tot if tot else 0
    out.write(f"{tag}\tintron_support\t{sup}/{tot}\t{pct:.2f}%\n")
PY
BASE="$OUTDIR/03_junctions"

# 4) Suporte de exons por cobertura (>0 em pelo menos um BAM)
for tag in noG stB stBH; do
  GTF="${MAP[$tag]}"
  awk 'BEGIN{OFS="\t"} $3=="exon"{print $1,$4-1,$5,$9,$6,$7}' "$GTF" > "$OUTDIR/04_exon_support/${tag}.exons.bed"
  bedtools multicov -bams "$BAM_DIR"/*.sorted.bam -bed "$OUTDIR/04_exon_support/${tag}.exons.bed" \
    > "$OUTDIR/04_exon_support/${tag}.exons.cov"
  awk '{c=0; for(i=7;i<=NF;i++) if($i>0) {c=1; break}; if(c==1) sup++} END{print sup+0}' \
    "$OUTDIR/04_exon_support/${tag}.exons.cov" > "$OUTDIR/04_exon_support/${tag}.exons_supported.count"
  wc -l "$OUTDIR/04_exon_support/${tag}.exons.bed" | awk '{print $1}' > "$OUTDIR/04_exon_support/${tag}.exons_total.count"
done

python - <<'PY'
import os
base=os.environ["BASE"]
tags=["noG","stB","stBH"]
with open(os.path.join(base,"summary.txt"),"w") as out:
  for tag in tags:
    sup=int(open(os.path.join(base,f"{tag}.exons_supported.count")).read().strip())
    tot=int(open(os.path.join(base,f"{tag}.exons_total.count")).read().strip())
    pct=100.0*sup/tot if tot else 0
    out.write(f"{tag}\texon_cov>0\t{sup}/{tot}\t{pct:.2f}%\n")
PY
BASE="$OUTDIR/04_exon_support"

echo "✅ Pré-PASA pronto. Veja:"
echo " - $OUTDIR/03_junctions/summary.txt"
echo " - $OUTDIR/04_exon_support/summary.txt"
echo " - (opcional) BUSCO em $OUTDIR/02_busco/"
