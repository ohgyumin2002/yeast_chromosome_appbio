#merging forward and reverse counts

wt_1 <- read.csv("../deseq/wt_r1_1_counts.csv", row.names = 1)
wt_2 <- read.csv("../deseq/wt_r1_2_counts.csv", row.names = 1)

View(wt_1)

count_files <- list.files(path = "~/Desktop/yeast_chr/deseq/", pattern = "*_counts.csv", full.names = TRUE)
count_files

counts <- merge(wt_1, wt_2, by = "row.names", all = T)
rownames(counts) <- counts$Row.names
counts <- counts[, -1]
View(counts)

# Loop through each counts file and load into R
for (file in count_files) {
  # Read the count file
  data <- read.csv(file, row.names = 1)
  
  counts <- merge(counts, data, by = "row.names", all = T)
  rownames(counts) <- counts$Row.names
  counts <- counts[, -1]
  View(counts)
}

counts <- counts[, c(-1, -2)]
View(counts)

all_counts <- counts[, c(7,8,9,10,11,12,1,2,3,4,5,6)]
View(all_counts)

write.csv(all_counts, file = "../deseq/allCounts.csv")