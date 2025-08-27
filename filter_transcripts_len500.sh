#!/usr/bin/env bash
set -euo pipefail

GFF="${1:-}"
THRESH="${2:-500}"
FASTA="${3:-}"   # opcional

if [[ -z "${GFF}" ]]; then
  echo "Uso: $0 <anotacoes.gff3> [limiar_bp=500] [genoma.fa(opcional)]" >&2
  exit 1
fi
if [[ ! -f "${GFF}" ]]; then
  echo "ERRO: GFF3 '${GFF}' não encontrado." >&2
  exit 1
fi

echo ">> Arquivo GFF3: ${GFF}"
echo ">> Limiar: ${THRESH} pb"
[[ -n "${FASTA}" ]] && echo ">> Genoma FASTA: ${FASTA}"

workdir="$(mktemp -d)"
trap 'rm -rf "${workdir}"' EXIT

############################################
# OPÇÃO B: se FASTA foi fornecido (gffread)
############################################
if [[ -n "${FASTA}" ]]; then
  if ! command -v gffread >/dev/null 2>&1; then
    echo "ERRO: 'gffread' não encontrado no PATH. Instale (conda: bioconda:gffread) ou rode sem FASTA." >&2
    exit 1
  fi
  if ! command -v seqkit >/dev/null 2>&1; then
    echo "ERRO: 'seqkit' não encontrado no PATH. Instale (bioconda:seqkit) ou rode sem FASTA." >&2
    exit 1
  fi

  echo ">> [B] Extraindo sequências de transcritos com gffread..."
  # -w escreve sequências spliced por transcrito; cabeçalho tende a ser o ID do transcrito
  gffread -w "${workdir}/transcripts.fa" -g "${FASTA}" "${GFF}"

  echo ">> [B] Medindo comprimentos e aplicando corte >= ${THRESH} pb..."
  # fx2tab: -n = nome, -l = length
  # garantimos que o ID seja apenas o primeiro token do defline
  seqkit fx2tab -n -l "${workdir}/transcripts.fa" \
    | awk -v T="${THRESH}" '{
        split($1,a," "); id=a[1]; len=$2;
        if (len>=T) print id;
      }' \
    | sort -u > "${workdir}/keep_transcripts.txt"

else
############################################
# OPÇÃO A: sem FASTA (somar exons / CDS)
############################################
  echo ">> [A] Calculando comprimentos por transcrito a partir de exons (fallback: CDS) ..."
  # Somar exons por Parent (transcrito); se não houver exon, usar soma de CDS
  awk -F'\t' 'BEGIN{OFS="\t"}
    $0 ~ /^#/ {next}
    {
      feat=$3; start=$4; end=$5; attr=$9;
      len = (end - start + 1);
      if (match(attr,/Parent=([^;]+)/,p)) {
        n=split(p[1], pp, /,/);
        if (feat=="exon") {
          for (i=1;i<=n;i++) exon_len[pp[i]] += len;
        } else if (feat=="CDS") {
          for (i=1;i<=n;i++) cds_len[pp[i]] += len;
        }
      }
    }
    END{
      # preferir exon; se não houver, usar cds
      for (t in exon_len) print t, exon_len[t] > "'"${workdir}"'/exon_lengths.tsv";
      for (t in cds_len)  print t, cds_len[t]  > "'"${workdir}"'/cds_lengths.tsv";
    }' "${GFF}"

  # construir tabela final transcript_id \t length
  if [[ -s "${workdir}/exon_lengths.tsv" ]]; then
    awk 'NR==FNR{ex[$1]=$2;next} {if(!($1 in ex)) ex[$1]=$2} END{for(i in ex) print i"\t"ex[i]}' \
      "${workdir}/exon_lengths.tsv" "${workdir}/cds_lengths.tsv" \
      | LC_ALL=C sort -k1,1 > "${workdir}/transcript_lengths.tsv"
  else
    # sem exons — cair para CDS apenas
    LC_ALL=C sort -k1,1 "${workdir}/cds_lengths.tsv" > "${workdir}/transcript_lengths.tsv"
  fi

  # manter apenas transcritos >= THRESH
  awk -v T="${THRESH}" '$2>=T{print $1}' "${workdir}/transcript_lengths.tsv" \
    | sort -u > "${workdir}/keep_transcripts.txt"
fi

# Se nada passou, sair cedo
if [[ ! -s "${workdir}/keep_transcripts.txt" ]]; then
  echo "⚠️  Nenhum transcrito >= ${THRESH} pb encontrado." >&2
  exit 0
fi

# Obter genes-parents dos transcritos mantidos
awk -F'\t' 'NR==FNR{keep[$1]=1; next}
  $0 ~ /^#/ {next}
  $3 ~ /(mRNA|transcript|lnc_RNA|ncRNA|rRNA|tRNA|snRNA|snoRNA)/ {
    id=""; par="";
    if (match($9,/ID=([^;]+)/,a)) id=a[1];
    if (match($9,/Parent=([^;]+)/,b)) par=b[1];
    if (keep[id] && par!="") {
      n=split(par,pp,/,/);
      for(i=1;i<=n;i++) print pp[i];
    }
  }' "${workdir}/keep_transcripts.txt" "${GFF}" \
  | sort -u > "${workdir}/keep_genes.txt"

# Conjunto de IDs a manter (genes + transcritos)
cat "${workdir}/keep_transcripts.txt" "${workdir}/keep_genes.txt" \
  | sort -u > "${workdir}/keep_ids.txt"

# Filtrar GFF: manter genes (ID em keep), transcritos (ID em keep) e
# quaisquer features filhas (Parent em keep de transcrito)
awk -F'\t' 'BEGIN{OFS="\t"}
  NR==FNR{keep[$1]=1; next}
  {
    if ($0 ~ /^#/) { print; next }
    feat=$3; attrs=$9; id=""; parent="";
    if (match(attrs,/ID=([^;]+)/,a)) id=a[1];
    if (match(attrs,/Parent=([^;]+)/,b)) parent=b[1];

    keep_line=0
    # manter gene/transcrito cujo ID está em keep
    if (feat=="gene" && id!="" && (id in keep)) keep_line=1
    if (feat ~ /(mRNA|transcript|lnc_RNA|ncRNA|rRNA|tRNA|snRNA|snoRNA)/ && id!="" && (id in keep)) keep_line=1

    # manter features cujos pais (Parent) incluem transcritos mantidos
    if (parent!=""){
      n=split(parent,pp,/,/);
      for(i=1;i<=n;i++){
        if (pp[i] in keep){ keep_line=1; break }
      }
    }

    if (keep_line) print
  }' "${workdir}/keep_ids.txt" "${GFF}" > filtered_len${THRESH}.gff3

# Pequeno sumário
kept_t=$(wc -l < "${workdir}/keep_transcripts.txt" | tr -d ' ')
kept_g=$(wc -l < "${workdir}/keep_genes.txt" | tr -d ' ')
echo "✅ Arquivo gerado: filtered_len${THRESH}.gff3"
echo "   Transcritos mantidos: ${kept_t}"
echo "   Genes correspondentes: ${kept_g}"
