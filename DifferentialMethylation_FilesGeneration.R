setwd("/home/zarifin_1/methylation")

#libraries
library(tidyverse)
library(GenomicRanges)
library(AnnotationHub)
library(IRanges)
library(ggrepel)
library(cowplot)
library(openxlsx)

#GTF file input for reference (GENCODE Release 23 (GRCh38.p3))
# ah <- AnnotationHub() # some problem with download, GTF reference was downloaded from local instead and put into data/GC_AH75128.rds
# gc <- ah[["AH75128"]] 
gc <- readRDS("data/GC_AH75128.rds")

# AH49556 hg38
# AH75128 hg19 (liftover)

#input files
soft <- read_table("data/soft_BS-dedup-BS_of_CpGDiad_bismark_bt2-CpG_All_BothStrands_Unique.txt", col_names = F)
stiff <- read_table("data/stiff_BS-dedup-BS_of_CpGDiad_bismark_bt2-CpG_All_BothStrands_Unique copy.txt", col_names = F)
colnames(soft) <- c("index", "chr", "loci", "strand", "context", "methylation")
colnames(stiff) <- c("index", "chr", "loci", "strand", "context", "methylation")

## problem: why does the chromosome goes from 1 to 27? And there is no X and/or Y.

# files prep
df_list <- list(soft = soft, stiff = stiff)

#change 23 to X and 24 to Y
for (i in 1:length(df_list)){
  df_list[[i]]$chr <- ifelse(df_list[[i]]$chr ==  23, "X", 
                             ifelse(df_list[[i]]$chr ==  24, "Y", 
                                    ifelse(df_list[[i]]$chr ==  25, "M",df_list[[i]]$chr
                                           )))
}

grange_list <- as.list(rep(NA, length(df_list)))
overlap_list <- as.list(rep(NA, length(df_list)))
genes <- as.list(rep(NA, length(df_list)))

for (i in 1:length(df_list)){
  df_list[[i]]$chr <- paste0("chr",df_list[[i]]$chr)
  df_list[[i]]$loci2 <- df_list[[i]]$loci + 0
  grange_list[[i]] <- df_list[[i]] %>%
    dplyr::select(chr, loci, loci2, strand) %>%
    dplyr::rename(chr = chr, start = loci, end = loci2, strand = strand)
  grange_list[[i]]$strand <- ifelse(grange_list[[i]]$strand %in% 1, "+", "-")
  # 1 is changed to "+", -1 is changed to "-"
  grange_list[[i]] <- makeGRangesFromDataFrame(grange_list[[i]])

  overlap_list[[i]] <- findOverlaps(grange_list[[i]], gc)
  genes[[i]] <- IRanges::extractList(gc$gene_name, as(overlap_list[[i]], "List"))
  genes[[i]] <- unstrsplit(unique(genes[[i]]), ";") # Needed in case more than one gene overlaps.

  df_list[[i]]$gene <- genes[[i]]
  
  grange_list[[i]]$gene <- genes[[i]]
}
names(grange_list) <- c("soft", "stiff")
saveRDS(grange_list, "results/grange_list.rds")
saveRDS(df_list, "results/df_list_with_gene.rds")

ratio_table <-  as.list(rep(NA, length(df_list)))

# aggregation
for (i in 1:length(df_list)){
  ratio_table[[i]] <- df_list[[i]] %>%
    group_by(chr, gene) %>%
    dplyr::count(methylation) %>%
    ungroup()
  
  ratio_table[[i]] <- ratio_table[[i]][!(is.na(ratio_table[[i]]$gene) | ratio_table[[i]]$gene==""), ] #remove no gene
  
  ratio_table[[i]] <- ratio_table[[i]] %>%  
    group_by(gene) %>% 
    summarise(prop = n[methylation == "Z"] / (n[methylation == "z"] + n[methylation == "Z"]), 
              n_gene =  n[methylation == "Z"] + n[methylation == "z"])
}
names(ratio_table) <- c("soft", "stiff")
openxlsx::write.xlsx(ratio_table, "results/ratio_table.xlsx")


