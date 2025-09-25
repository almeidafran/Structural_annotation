#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
[03] Final EVM: QC completo

Uso:
  03_finalEVM_qc.sh \
    -g <genome.fa> \
    -a <evm_final.gff3> \
    -b <bam_dir> \
    -t <threads> \
    -l <busco_lineage> \
    -o <outdir>

Saídas:
  <outdir>/
    01_stats/agat_evm.txt
    02_busco/ (transcriptome + proteins)
    03_junctions/summary.txt
    04_exon_support/summary.txt
USAGE
}

GENOME="" ; EVM_GFF="" ; BAM_DIR="" ; THREADS=8 ; LINEAGE="" ; OUTDIR=""
while getopts ":g:a:b:t:l:o:h" opt; do
  case "$opt" in
    g) GENOME="$OPTARG";;
    a) EVM_GFF="$OPTARG";;
    b) BAM_DIR="$OPTARG";;
    t) THREADS="$OPTARG";;
    l) LINEAGE="$OPTARG";;
    o) OUTDIR="$OPTARG";;
    h) usage; exit 0;;
  esac
done
[[ -f "$GENOME" && -f "$EVM_GFF" && -d "$BAM_DIR" && -n "$OUTDIR" && -n "$LINEAGE" ]] || { usage; exit 2; }

for bin in agat_sp_statistics.pl gffread busco bedtools regtools; do
  command -v "$bin" >/dev/null 2>&1 || { echo "❌ falta $bin"; exit 1; }
done

mkdir -p "$OUTDIR"/{01_stats,02_busco,03_junctions/junctions_bed12,04_exon_support,logs}

# 1) AGAT
agat_sp_statistics.pl --gff "$EVM_GFF" -g "$GENOME" -o "$OUTDIR/01_stats/agat_evm.txt" > "$OUTDIR/logs/agat_evm.log" 2>&1

# 2) BUSCO em transcritos e proteínas
gffread "$EVM_GFF" -g "$GENOME" \
  -w "$OUTDIR/02_busco/evm.transcripts.fa" \
  -x "$OUTDIR/02_busco/evm.cds.fa" \
  -y "$OUTDIR/02_busco/evm.proteins.fa"
busco -i "$OUTDIR/02_busco/evm.transcripts.fa" -l "$LINEAGE" -m transcriptome -c "$THREADS" \
  -o busco_tx_evm --out_path "$OUTDIR/02_busco" > "$OUTDIR/logs/busco_tx_evm.log" 2>&1 || true
busco -i "$OUTDIR/02_busco/evm.proteins.fa" -l "$LINEAGE" -m proteins -c "$THREADS" \
  -o busco_prot_evm --out_path "$OUTDIR/02_busco" > "$OUTDIR/logs/busco_prot_evm.log" 2>&1 || true

# 3) junções suportadas
for BAM in "$BAM_DIR"/*.sorted.bam; do
  S=$(basename "$BAM" .sorted.bam)
  regtools junctions extract -a 8 -m 50 -M 500000 "$BAM" > "$OUTDIR/03_junctions/junctions_bed12/${S}.bed12"
done
cat "$OUTDIR"/03_junctions/junctions_bed12/*.bed12 > "$OUTDIR/03_junctions/all_junctions.bed12"
agat_convert_sp_gff2bed.pl --gff "$EVM_GFF" --out "$OUTDIR/03_junctions/evm.bed" > /dev/null 2>&1
awk '$8=="intron"{OFS="\t"; print $1,$2,$3,$4,$5,$6}' "$OUTDIR/03_junctions/evm.bed" > "$OUTDIR/03_junctions/evm.introns.bed"
bedtools intersect -u -a "$OUTDIR/03_junctions/evm.introns.bed" -b "$OUTDIR/03_junctions/all_junctions.bed12" \
  | wc -l > "$OUTDIR/03_junctions/introns_supported.count"
wc -l "$OUTDIR/03_junctions/evm.introns.bed" | awk '{print $1}' > "$OUTDIR/03_junctions/introns_total.count"
python - <<'PY'
import os
b=os.environ["B"]
sup=int(open(os.path.join(b,"introns_supported.count")).read().strip())
tot=int(open(os.path.join(b,"introns_total.count")).read().strip())
pct=100.0*sup/tot if tot else 0
open(os.path.join(b,"summary.txt"),"w").write(f"evm\tintron_support\t{sup}/{tot}\t{pct:.2f}%\n")
PY
B="$OUTDIR/03_junctions"

# 4) exons com cobertura > 0
awk 'BEGIN{OFS="\t"} $3=="exon"{print $1,$4-1,$5,$9,$6,$7}' "$EVM_GFF" > "$OUTDIR/04_exon_support/evm.exons.bed"
bedtools multicov -bams "$BAM_DIR"/*.sorted.bam -bed "$OUTDIR/04_exon_support/evm.exons.bed" \
  > "$OUTDIR/04_exon_support/evm.exons.cov"
awk '{c=0; for(i=7;i<=NF;i++) if($i>0){c=1;break}; if(c==1) sup++} END{print sup+0}' \
  "$OUTDIR/04_exon_support/evm.exons.cov" > "$OUTDIR/04_exon_support/evm.exons_supported.count"
wc -l "$OUTDIR/04_exon_support/evm.exons.bed" | awk '{print $1}' > "$OUTDIR/04_exon_support/evm.exons_total.count"
python - <<'PY'
import os
b=os.environ["B"]
sup=int(open(os.path.join(b,"evm.exons_supported.count")).read().strip())
tot=int(open(os.path.join(b,"evm.exons_total.count")).read().strip())
pct=100.0*sup/tot if tot else 0
open(os.path.join(b,"summary.txt"),"w").write(f"evm\texon_cov>0\t{sup}/{tot}\t{pct:.2f}%\n")
PY
B="$OUTDIR/04_exon_support"

echo "✅ QC final do EVM pronto em: $OUTDIR"
