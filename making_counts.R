setwd("~/Desktop/yeast_chr/bam_files/")

library(limma)

library(Rsubread)
fc <- featureCounts(files = "sy14_r3_2.sorted.bam", 
                    annot.ext = "~/Desktop/ncbi_dataset/ncbi_dataset/data/GCF_000146045.2/genomic.gtf", 
                    isGTFAnnotationFile = TRUE,
                    isPairedEnd = FALSE)
View(fc$counts)
write.csv(fc$counts, "~/Desktop/yeast_chr/deseq/sy14_r3_2_counts.csv")
