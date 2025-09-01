#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# make_crl.sh
#  - Cria uma Custom Repeat Library (CRL) a partir da saída do RepeatModeler,
#    removendo consensos com similaridade proteica usando ProtExcluder.
#
# Requisitos:
#   - diamond (para blastx rápido)
#   - perl
#   - ProtExcluder.pl (indique via -x ou $PROTEX_PATH, tentamos ~/ProtExcluder/ProtExcluder.pl também)
#
# Uso:
#   ./make_crl.sh -r <repeats.fa> -p <proteins.faa> [-o <outprefix>] [-t <threads>] [-x <ProtExcluder.pl>]
#
# Exemplo:
#   ./make_crl.sh \
#       -r Vper.hapA.v1.unmasked_db-families.fa \
#       -p /home/externo/francined/proteomes_swiss_prot/swissprot_plants.min50.faa \
#       -o Vper.hapA.CRL_from_ProtExcluder \
#       -t 24
#
# Saídas principais:
#   <outprefix>/
#     ├─ diamond_db/                             (db do diamond)
#     ├─ repeats_vs_proteins.tab                 (resultados do blastx)
#     ├─ CRL.filtered.fa                         (biblioteca filtrada = sua CRL)
#     ├─ kept_ids.txt                            (IDs mantidos)
#     ├─ removed_by_ProtExcluder.txt             (IDs removidos)
#     └─ summary.txt                             (resumo Antes/Depois)
# ------------------------------------------------------------

usage() {
  cat <<EOF
Uso: $(basename "$0") -r <repeats.fa> -p <proteins.faa> [-o <outprefix>] [-t <threads>] [-x <ProtExcluder.pl>]

  -r  FASTA de repetições (RepeatModeler), ex: consensi.fa.classified
  -p  FASTA de proteínas de referência (Swiss-Prot plantas, etc.)
  -o  Prefixo/dir de saída (default: basename do repeats + ".CRL")
  -t  Threads (default: 8)
  -x  Caminho do ProtExcluder.pl (opcional)
  -h  Esta ajuda
EOF
}

REPEATS=""
PROTEINS=""
OUTPREFIX=""
THREADS=8
PROTEX_PATH="${PROTEX_PATH:-}"

while getopts ":r:p:o:t:x:h" opt; do
  case $opt in
    r) REPEATS="$OPTARG" ;;
    p) PROTEINS="$OPTARG" ;;
    o) OUTPREFIX="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    x) PROTEX_PATH="$OPTARG" ;;
    h) usage; exit 0 ;;
    \?) echo "Opção inválida: -$OPTARG" >&2; usage; exit 2 ;;
    :) echo "A opção -$OPTARG requer um argumento." >&2; usage; exit 2 ;;
  case
  esac
done

# Checagens básicas
if [[ -z "${REPEATS}" || -z "${PROTEINS}" ]]; then
  echo "ERRO: informe -r <repeats.fa> e -p <proteins.faa>." >&2
  usage
  exit 2
fi
if [[ ! -s "$REPEATS" ]]; then
  echo "ERRO: arquivo de repetições inexistente/vazio: $REPEATS" >&2
  exit 1
fi
if [[ ! -s "$PROTEINS" ]]; then
  echo "ERRO: arquivo de proteínas inexistente/vazio: $PROTEINS" >&2
  exit 1
fi

# Ferramentas
command -v diamond >/dev/null 2>&1 || { echo "ERRO: 'diamond' não encontrado no PATH." >&2; exit 1; }
command -v perl    >/dev/null 2>&1 || { echo "ERRO: 'perl' não encontrado no PATH."    >&2; exit 1; }

# Localizar ProtExcluder.pl
if [[ -z "${PROTEX_PATH}" ]]; then
  if [[ -x "$HOME/ProtExcluder/ProtExcluder.pl" ]]; then
    PROTEX_PATH="$HOME/ProtExcluder/ProtExcluder.pl"
  elif command -v ProtExcluder.pl >/dev/null 2>&1; then
    PROTEX_PATH="$(command -v ProtExcluder.pl)"
  else
    echo "ERRO: ProtExcluder.pl não encontrado. Informe via -x <caminho> ou PROTEX_PATH." >&2
    exit 1
  fi
fi
if [[ ! -s "$PROTEX_PATH" ]]; then
  echo "ERRO: ProtExcluder.pl inválido: $PROTEX_PATH" >&2
  exit 1
fi

# Prefixo de saída
if [[ -z "${OUTPREFIX}" ]]; then
  base="$(basename "$REPEATS")"
  OUTPREFIX="${base%.*}.CRL"
fi
OUTDIR="$OUTPREFIX"
mkdir -p "$OUTDIR"

echo ">>> Repeats:     $REPEATS"
echo ">>> Proteins:    $PROTEINS"
echo ">>> ProtExcluder: $PROTEX_PATH"
echo ">>> Threads:     $THREADS"
echo ">>> Outdir:      $OUTDIR"
echo

# 1) Preparar DB DIAMOND
DBDIR="$OUTDIR/diamond_db"
mkdir -p "$DBDIR"
DB="$DBDIR/proteins"
if [[ ! -s "${DB}.dmnd" ]]; then
  echo "[1/4] Criando DB DIAMOND..."
  diamond makedb -p "$THREADS" -in "$PROTEINS" -d "$DB"
else
  echo "[1/4] DB DIAMOND já existe, pulando ( ${DB}.dmnd )"
fi

# 2) DIAMOND blastx (repeats vs proteins)
BLAST_TAB="$OUTDIR/repeats_vs_proteins.tab"
echo "[2/4] Rodando DIAMOND blastx..."
diamond blastx \
  -d "$DB" \
  -q "$REPEATS" \
  -o "$BLAST_TAB" \
  -e 1e-5 \
  -k 1 \
  --outfmt 6 \
  --threads "$THREADS"

# 3) Rodar ProtExcluder
echo "[3/4] Rodando ProtExcluder..."
# ProtExcluder costuma escrever a(s) saídas no diretório de trabalho atual;
# vamos rodar dentro do OUTDIR para manter tudo organizado.
(
  cd "$OUTDIR"
  # usamos caminhos relativos dentro do outdir
  perl "$PROTEX_PATH" "$(basename "$BLAST_TAB")" "$(realpath "$REPEATS")"
)

# 4) Detectar arquivo filtrado e gerar sumário
echo "[4/4] Gerando resumo e nomes padronizados..."

# Contagem original
ORIG_N=$(grep -c '^>' "$REPEATS" || true)

# Tentar localizar o FASTA filtrado gerado pelo ProtExcluder
# Procuramos por arquivos novos dentro do OUTDIR que pareçam FASTA e não sejam o original.
FILTERED_CANDIDATES=$(ls -1 "$OUTDIR"/* 2>/dev/null | grep -Ei '\.(fa|fasta|fas)(\.|$)' || true)
FILTERED=""
for f in $FILTERED_CANDIDATES; do
  # ignore o BLAST tab e outros
  [[ "$f" == "$BLAST_TAB" ]] && continue
  # se for idêntico ao repeats original (mesmo arquivo), ignore
  if cmp -s "$f" "$REPEATS"; then
    continue
  fi
  # heurística: o filtrado deve ter <= ao número de sequências original
  n=$(grep -c '^>' "$f" || true)
  if [[ -n "$n" && "$n" -le "$ORIG_N" && "$n" -gt 0 ]]; then
    FILTERED="$f"
    break
  fi
done

if [[ -z "$FILTERED" ]]; then
  echo "ATENÇÃO: não consegui identificar automaticamente o FASTA filtrado do ProtExcluder no diretório $OUTDIR."
  echo "Verifique os arquivos gerados e renomeie manualmente o filtrado para '$OUTDIR/CRL.filtered.fa'."
  exit 0
fi

# Padronizar nome do filtrado
CRL="$OUTDIR/CRL.filtered.fa"
cp -f "$FILTERED" "$CRL"

# IDs mantidos/removidos
grep -E '^>' "$CRL"   | sed 's/^>//' | sort > "$OUTDIR/kept_ids.txt"
grep -E '^>' "$REPEATS" | sed 's/^>//' | sort > "$OUTDIR/all_ids.txt"
comm -23 "$OUTDIR/all_ids.txt" "$OUTDIR/kept_ids.txt" > "$OUTDIR/removed_by_ProtExcluder.txt" || true

KEPT_N=$(wc -l < "$OUTDIR/kept_ids.txt" | tr -d '[:space:]')
REMV_N=$(wc -l < "$OUTDIR/removed_by_ProtExcluder.txt" | tr -d '[:space:]')

# Sumário
cat > "$OUTDIR/summary.txt" <<EOS
Custom Repeat Library (ProtExcluder)
------------------------------------
Repetições de entrada  : $REPEATS
Proteínas (referência) : $PROTEINS
BLAST tab              : $BLAST_TAB
CRL (filtrada)         : $CRL

Contagem:
  Antes (consensos): $ORIG_N
  Mantidos         : $KEPT_N
  Removidos        : $REMV_N
EOS

echo
echo "✅ Concluído."
echo "  - CRL:                $CRL"
echo "  - BLAST table:        $BLAST_TAB"
echo "  - kept_ids.txt:       $OUTDIR/kept_ids.txt"
echo "  - removed_by_ProtExcluder.txt: $OUTDIR/removed_by_ProtExcluder.txt"
echo "  - summary.txt:        $OUTDIR/summary.txt"
