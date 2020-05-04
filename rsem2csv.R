# Load modules
library(Tmisc)
library(annotables)
library(tidyverse)
library(tximport)

source("./collapseIsoforms.R")
source("./fixERCCs.R")

study <- "tutorial"

# Collapse the annotation table to the ensembl gene id (mouse)
annotables::grcm38 %>% distinct(ensgene)
anno_collapsed_to_ensgene <- annotables::grcm38 %>% 
  group_by(ensgene) %>% 
  summarize_all(funs(. %>% unique %>% paste(collapse="; ")))

sampnames <- list.files(path = "../3-rsem",recursive = TRUE,pattern = "isoforms.results") %>% str_split_fixed("/",n = 2)
sampnames <- sampnames[,1]
quantfiles <- list.files(path = "../3-rsem",recursive = TRUE,pattern = "isoforms.results", full.names = T) %>% set_names(sampnames)

# Join mouse genes and ERCCs
tx2gene <- annotables::grcm38_tx2gene
ercctx2gene <- read_tsv("./ercc_tx2gene.tsv")
tx2gene <- rbind(tx2gene, ercctx2gene)

# Import transcript data
txi <- tximport(quantfiles,
                type = "rsem",
                tx2gene = tx2gene,
                importer = read_tsv,
                txIn = TRUE,
                txOut = FALSE,
                ignoreTxVersion = TRUE, 
                countsFromAbundance="no")

# Annotate genes
rsem.counts <- txi$counts %>% as.data.frame()
rsem.counts$ensgene <- rownames(rsem.counts)
rsem.counts <- left_join(rsem.counts, anno_collapsed_to_ensgene) %>%
  select(ensgene, entrez, symbol, chr, start, end, strand, biotype, description, everything())
index_counts <- 10:ncol(rsem.counts)

# Fix ERCC annotation
rsem.counts <- fixERCCs(rsem.counts)

# Collapse isoforms to genes
rsem.counts <- collapseIsoforms(rsem.counts,index_counts)

rsem.counts %>% write_csv(paste0("rsem_",study,".csv"))