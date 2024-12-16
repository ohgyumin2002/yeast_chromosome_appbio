res <- read.csv("deseq2_results.csv")
View(res)

deg <- res[res$padj<0.001 & (res$log2FoldChange >= 1 | res$log2FoldChange <= -1) , ]
deg <- na.omit(deg)
dim(deg)

gene_names <- deg$X
gene_names
write.csv(gene_names, "degs.csv")

degs <- read.csv("degs.csv")
View(degs)

length(gene_names)

degs_trans <- read.csv("translated_genes.csv")
dim(degs_trans)

deg$name <- degs_trans$name
View(deg)

upreg <- deg[deg$log2FoldChange < 0, ]
downreg <- deg[deg$log2FoldChange > 0, ]
dim(upreg)
dim(downreg)

compare_genes <- deg[deg$name %in% c("ERR2", "HSP32", "FEX2", "YPL277C", "MPH3",
                                     "VTH1", "SEO1", "YOL162W","YOL163W", "THI5",
                                     "MAL11", "YFR057W", "HSP26", "RNR3", "NCA3",
                                     "HBN1", "OYE3", "HUG1", "HSP12", "YKU80",
                                     "RDN37-1", "FIT3", "YBR012W-A", "YOL160W",
                                     "YDR261W-A", "TAR1", "YIR042C", "MET14"),]

dim(compare_genes)

View(compare_genes)