#!/usr/bin/env bash
set -euo pipefail

GUIDE="helixer_annotation.gff3"
THREADS=12
PACBIO_BAM="pacbio_isoseq.sorted.bam"
SR_BAMS=( VP1C.sorted.bam VP1F.sorted.bam VP1R.sorted.bam VP2C.sorted.bam VP2F.sorted.bam VP2R.sorted.bam VP3C.sorted.bam VP3F.sorted.bam VP3R.sorted.bam )

OUT=stringtie_per_sample
mkdir -p "$OUT/gtf_per_sample"

# 1) Montar cada SR BAM com guia
for B in "${SR_BAMS[@]}"; do
  S=$(basename "$B" .sorted.bam)
  stringtie -p "$THREADS" -G "$GUIDE" -o "$OUT/gtf_per_sample/${S}.gtf" "$B"
done

# 2) Montar LR (PacBio) com -L (long reads) + guia
stringtie -L -p "$THREADS" -G "$GUIDE" -o "$OUT/gtf_per_sample/PacBio.gtf" "$PACBIO_BAM"

# 3) Mesclar todas as montagens
ls "$OUT/gtf_per_sample/"*.gtf > "$OUT/mergelist.txt"
stringtie --merge -p "$THREADS" -G "$GUIDE" -o "$OUT/merged.gtf" "$OUT/mergelist.txt"

# 4) Avaliar contra a referência
gffcompare -r "$GUIDE" -o "$OUT/gffcmp" "$OUT/merged.gtf"

# 5) Quantificar por amostra usando a montagem mesclada
mkdir -p "$OUT/quant"
for B in "${SR_BAMS[@]}" "$PACBIO_BAM"; do
  S=$(basename "$B" .sorted.bam)
  stringtie -p "$THREADS" -G "$OUT/merged.gtf" -e -B \
    -A "$OUT/quant/${S}.abundance.tsv" \
    -o "$OUT/quant/${S}.gtf" "$B"
done
