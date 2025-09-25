#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
[02] Pós-PASA: QC do GFF/GTF refinado pelo PASA

Uso:
  02_posPASA_qc.sh \
    -g <genome.fa> \
    -a <pasa_refined.gff3> \
    -x <stringtie_escolhido.gtf> \
    -b <bam_dir> \
    -t <threads> \
    -o <outdir>

Saídas:
  <outdir>/
    01_stats/agat_pasa.txt
    02_compare/stringtie_vs_pasa.*  (gffcompare)
    03_junctions/summary.txt
    04_exon_support/summary.txt
USAGE
}

GENOME="" ; PASA_GFF="" ; ST_CHOSEN_GTF="" ; BAM_DIR="" ; THREADS=8 ; OUTDIR=""
while getopts ":g:a:x:b:t:o:h" opt; do
  case "$opt" in
    g) GENOME="$OPTARG";;
    a) PASA_GFF="$OPTARG";;
    x) ST_CHOSEN_GTF="$OPTARG";;
    b) BAM_DIR="$OPTARG";;
    t) THREADS="$OPTARG";;
    o) OUTDIR="$OPTARG";;
    h) usage; exit 0;;
  esac
done
[[ -f "$GENOME" && -f "$PASA_GFF" && -f "$ST_CHOSEN_GTF" && -d "$BAM_DIR" && -n "$OUTDIR" ]] || { usage; exit 2; }

for bin in agat_sp_statistics.pl gffcompare bedtools regtools; do
  command -v "$bin" >/dev/null 2>&1 || { echo "❌ falta $bin"; exit 1; }
done

mkdir -p "$OUTDIR"/{01_stats,02_compare,03_junctions/junctions_bed12,04_exon_support,logs}

# 1) AGAT no PASA
agat_sp_statistics.pl --gff "$PASA_GFF" -g "$GENOME" -o "$OUTDIR/01_stats/agat_pasa.txt" > "$OUTDIR/logs/agat_pasa.log" 2>&1

# 2) gffcompare (StringTie escolhido vs PASA)
gffcompare -r "$PASA_GFF" -o "$OUTDIR/02_compare/stringtie_vs_pasa" "$ST_CHOSEN_GTF" > "$OUTDIR/logs/gffcompare.log" 2>&1

# 3) junções suportadas (como no script 01)
for BAM in "$BAM_DIR"/*.sorted.bam; do
  S=$(basename "$BAM" .sorted.bam)
  regtools junctions extract -a 8 -m 50 -M 500000 "$BAM" > "$OUTDIR/03_junctions/junctions_bed12/${S}.bed12"
done
cat "$OUTDIR"/03_junctions/junctions_bed12/*.bed12 > "$OUTDIR/03_junctions/all_junctions.bed12"

agat_convert_sp_gff2bed.pl --gff "$PASA_GFF" --out "$OUTDIR/03_junctions/pasa.bed" > /dev/null 2>&1
awk '$8=="intron"{OFS="\t"; print $1,$2,$3,$4,$5,$6}' "$OUTDIR/03_junctions/pasa.bed" > "$OUTDIR/03_junctions/pasa.introns.bed"

bedtools intersect -u -a "$OUTDIR/03_junctions/pasa.introns.bed" -b "$OUTDIR/03_junctions/all_junctions.bed12" \
  | wc -l > "$OUTDIR/03_junctions/introns_supported.count"
wc -l "$OUTDIR/03_junctions/pasa.introns.bed" | awk '{print $1}' > "$OUTDIR/03_junctions/introns_total.count"

python - <<'PY'
import os
base=os.environ["BASE"]
sup=int(open(os.path.join(base,"introns_supported.count")).read().strip())
tot=int(open(os.path.join(base,"introns_total.count")).read().strip())
pct=100.0*sup/tot if tot else 0
open(os.path.join(base,"summary.txt"),"w").write(f"pasa\tintron_support\t{sup}/{tot}\t{pct:.2f}%\n")
PY
BASE="$OUTDIR/03_junctions"

# 4) exons com cobertura > 0
awk 'BEGIN{OFS="\t"} $3=="exon"{print $1,$4-1,$5,$9,$6,$7}' "$PASA_GFF" > "$OUTDIR/04_exon_support/pasa.exons.bed"
bedtools multicov -bams "$BAM_DIR"/*.sorted.bam -bed "$OUTDIR/04_exon_support/pasa.exons.bed" \
  > "$OUTDIR/04_exon_support/pasa.exons.cov"
awk '{c=0; for(i=7;i<=NF;i++) if($i>0){c=1;break}; if(c==1) sup++} END{print sup+0}' \
  "$OUTDIR/04_exon_support/pasa.exons.cov" > "$OUTDIR/04_exon_support/pasa.exons_supported.count"
wc -l "$OUTDIR/04_exon_support/pasa.exons.bed" | awk '{print $1}' > "$OUTDIR/04_exon_support/pasa.exons_total.count"

python - <<'PY'
import os
base=os.environ["BASE"]
sup=int(open(os.path.join(base,"pasa.exons_supported.count")).read().strip())
tot=int(open(os.path.join(base,"pasa.exons_total.count")).read().strip())
pct=100.0*sup/tot if tot else 0
open(os.path.join(base,"summary.txt"),"w").write(f"pasa\texon_cov>0\t{sup}/{tot}\t{pct:.2f}%\n")
PY
BASE="$OUTDIR/04_exon_support"

echo "✅ Pós-PASA QC pronto em: $OUTDIR"
