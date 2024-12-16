library(dplyr)

# Define the count files (which files to use)
count_files <- list.files(path = "~/Desktop/yeast_chr/deseq/", pattern = "*_counts.csv", full.names = TRUE)

# Initialise an empty list to store data frames
count_list <- list()
View(count_files)
# Loop through each counts file and load into R
for (file in count_files) {
  # Read the count file
  count_data <- read.table(file, header = TRUE, row.names = 1, sep = ",")
  
  # Extract the sample name from the file name
  sample_name <- sub("_counts.csv", "", basename(file))
  
  # Rename the counts column to the sample name
  colnames(count_data) <- c(sample_name)
  
  # Add the data frame to the list
  count_list[[sample_name]] <- count_data
}


# Combine all data frames by row names (gene IDs)
combined_counts <- Reduce(function(x, y) merge(x, y, by = "row.names", all = TRUE), count_list)

# Fix row names after merging
rownames(combined_counts) <- combined_counts$Row.names
combined_counts
combined_counts <- combined_counts[ , -1]

# Replace missing data with zeros (if any genes are not present in some samples)
#combined_counts[is.na(combined_counts)] <- 0

# Write the combined matrix to a csv file
write.csv(combined_counts, file = "~/Desktop/yeast_chr/deseq/round1.csv", quote = FALSE, row.names = TRUE)

setwd("~/Desktop/yeast_chr/deseq/")
round_x <- read.csv("round1.csv")
View(round_x)


