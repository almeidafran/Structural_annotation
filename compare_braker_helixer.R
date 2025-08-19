#!/usr/bin/env Rscript

suppressMessages({
  library(Biostrings)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(reshape2)
})

# --------- caminhos (ajuste se necessário) ----------
busco_braker  <- "braker/busco/short_summary.txt"
busco_helixer <- "helixer/busco/short_summary.txt"
prot_braker   <- "braker/proteins.fa"
prot_helixer  <- "helixer/proteins.fa"
agat_braker   <- "braker/braker.stats.txt"
agat_helixer  <- "helixer/helixer.stats.txt"

# --------- helpers ----------
read_busco <- function(file, label){
  # procura linha com "C:" (BUSCO v5)
  lines <- readLines(file, warn = FALSE)
  L <- lines[grep("C:", lines)][1]
  stopifnot(length(L) == 1)
  # extrai números na ordem: C S D F M (percentuais)
  # Ex: "C:95.6%[S:92.3%,D:3.3%],F:1.0%,M:3.4%,n:1614"
  num <- as.numeric(gsub("[^0-9.]", "", unlist(strsplit(L, "[,\\[\\]]"))))
  # mapeia: C, S, D, F, M (desconsidera 'n')
  df <- data.frame(Method = label,
                   Complete   = num[1],
                   Single     = num[2],
                   Duplicated = num[3],
                   Fragmented = num[4],
                   Missing    = num[5],
                   check.names = FALSE)
  df
}

read_fasta_lengths <- function(file){
  if(!file.exists(file)) stop(paste("FASTA não encontrado:", file))
  w <- width(readAAStringSet(file))
  as.integer(w)
}

# tenta extrair % monoéxons e média de éxons do arquivo de stats do AGAT
parse_agat_stats <- function(file, label){
  if(!file.exists(file)) stop(paste("AGAT stats não encontrado:", file))
  lines <- readLines(file, warn = FALSE)

  # heurísticas para achar números (AGAT costuma imprimir em texto)
  mono_line <- lines[grepl("monoexon", tolower(lines))][1]
  exons_line <- lines[grepl("exon per (gene|mrna)|mean exon", tolower(lines))][1]

  mono <- NA_real_; exons <- NA_real_
  if(length(mono_line) == 1) {
    nums <- as.numeric(gsub(",", ".", regmatches(mono_line, gregexpr("[0-9]+\\.?[0-9]*", mono_line))[[1]]))
    mono <- nums[length(nums)]
  }
  if(length(exons_line) == 1) {
    nums <- as.numeric(gsub(",", ".", regmatches(exons_line, gregexpr("[0-9]+\\.?[0-9]*", exons_line))[[1]]))
    exons <- nums[length(nums)]
  }
  data.frame(Method = label,
             monoexon_percentage = mono,
             mean_exons_per_gene = exons)
}

# --------- BUSCO ----------
busco <- rbind(
  read_busco(busco_braker,  "BRAKER"),
  read_busco(busco_helixer, "Helixer")
)

busco_long <- melt(busco, id.vars = "Method")
p1 <- ggplot(busco_long, aes(x = Method, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Percentual (%)", x = NULL, fill = NULL, title = "BUSCO (proteínas)") +
  theme_bw(base_size = 12)
ggsave("busco_comparison.png", p1, width = 7, height = 4, dpi = 300)

# --------- Comprimento de proteínas ----------
lb <- read_fasta_lengths(prot_braker)
lh <- read_fasta_lengths(prot_helixer)
prot_df <- data.frame(Length = c(lb, lh),
                      Method = rep(c("BRAKER","Helixer"), c(length(lb), length(lh))))
p2 <- ggplot(prot_df, aes(x = Method, y = Length, fill = Method)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.4) +
  scale_y_continuous(trans = "log10") +
  labs(y = "Tamanho da proteína (aa, escala log10)", x = NULL, title = "Distribuição de tamanhos das proteínas") +
  theme_bw(base_size = 12)
ggsave("protein_length_comparison.png", p2, width = 7, height = 4, dpi = 300)

# --------- Monoéxons e éxons/gene (AGAT) ----------
agat <- rbind(
  parse_agat_stats(agat_braker,  "BRAKER"),
  parse_agat_stats(agat_helixer, "Helixer")
)

# Monoéxons
p3 <- ggplot(agat, aes(x = Method, y = monoexon_percentage, fill = Method)) +
  geom_bar(stat = "identity") +
  labs(y = "% genes monoéxons", x = NULL, title = "Genes monoéxons") +
  theme_bw(base_size = 12)
ggsave("monoexon_comparison.png", p3, width = 5.5, height = 4, dpi = 300)

# Éxons por gene
p4 <- ggplot(agat, aes(x = Method, y = mean_exons_per_gene, fill = Method)) +
  geom_bar(stat = "identity") +
  labs(y = "Éxons por gene (média)", x = NULL, title = "Estrutura gênica") +
  theme_bw(base_size = 12)
ggsave("exons_per_gene_comparison.png", p4, width = 5.5, height = 4, dpi = 300)

cat("✅ Gráficos gerados:\n",
    " - busco_comparison.png\n",
    " - protein_length_comparison.png\n",
    " - monoexon_comparison.png\n",
    " - exons_per_gene_comparison.png\n")
