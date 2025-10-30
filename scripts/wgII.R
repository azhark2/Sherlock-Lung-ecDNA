# This implementation:
#   1. Takes a seg_data dataframe with columns: chromosome, start, end, cn
# 2. Iterates through each chromosome
# 3. For each chromosome:
#   - Sets negative copy numbers to 0
# - Calculates ploidy as length-weighted average of copy numbers
# - Identifies segments with gains (CN > ploidy) and losses (CN < ploidy)
# - Calculates proportion of chromosome that is aberrant
# 4. Sums up the proportions and divides by 22 to get final wGII score

#input CNVkit calls
setwd("~/iCloud/dev/Sherlock-Lung/data/ShatterSeek-input/")
CN_files <- list.files(path = "~/iCloud/dev/Sherlock-Lung/data/ShatterSeek-input/CNVkit/", 
                       pattern = ".segments.tsv", 
                       full.names = FALSE)

calculate_wGII <- function(seg_data, chromosomes = paste0('chr', 1:22)) {
  wGII <- 0
  
  # Rename columns to match function expectations
  names(seg_data) <- c("chromosome", "start", "end", "cn")
  
  # Add chr prefix if not present
  seg_data$chromosome <- paste0("chr", seg_data$chromosome)
  
  # Set negative CN values to 0
  seg_data$cn[seg_data$cn < 0] <- 0
  
  for(chrom in chromosomes) {
    # Calculate ploidy for current data
    ploidy <- with(seg_data, sum((end - start)*cn) / sum(end - start))
    
    # Get segments for current chromosome
    sub_chr <- seg_data[seg_data$chromosome == chrom,]
    
    if(nrow(sub_chr) > 0) {  # Only process if chromosome has data
      # Calculate gains and losses
      gains <- round(sub_chr$cn) > round(ploidy)
      losses <- round(sub_chr$cn) < round(ploidy)
      
      seg_gain <- sum(abs(sub_chr$end[gains] - sub_chr$start[gains]))
      seg_loss <- sum(abs(sub_chr$end[losses] - sub_chr$start[losses]))
      seg_total <- sum(abs(sub_chr$end - sub_chr$start))
      
      # Calculate GII for this chromosome
      if(seg_total > 0) {  # Avoid division by zero
        GII <- (seg_gain + seg_loss) / seg_total
        wGII <- wGII + GII
      }
    }
  }
  
  # Average over all chromosomes
  wGII <- wGII / 22
  
  return(wGII)
}

# Process all files and store results
results_list <- list()
for (file in CN_files) {
  # Read CNVkit file
  cnv_data <- read.table(paste0("~/iCloud/dev/Sherlock-Lung/data/ShatterSeek-input/CNVkit/", file), 
                         header=TRUE, 
                         sep="\t",
                         stringsAsFactors=FALSE)
  
  sample <- gsub(".segments.tsv", "", file)
  print(paste("Processing sample:", sample))
  
  # Calculate wGII
  result <- tryCatch({
    calculate_wGII(cnv_data)
  }, error = function(e) {
    print(paste("Error processing", sample, ":", e$message))
    return(NA)
  })
  
  results_list[[sample]] <- result
}

# Convert results to data frame
results_df <- data.frame(
  sample = names(results_list),
  wGII = unlist(results_list),
  stringsAsFactors = FALSE
)

# Print results
print(results_df)

#rename sample column to Tumor_Barcode
colnames(results_df)[1] <- "Tumor_Barcode"

#write out to text file
write.table(results_df, "~/iCloud/dev/Sherlock-Lung/results2/wGII_results.txt", row.names = FALSE, quote = FALSE, sep = "\t")

#merge with log-reg4
data <- read.table("~/iCloud/dev/Sherlock-Lung/results2/log-reg-table4.tsv", header = TRUE, sep = "\t")
merged <- merge(data, results_df, by = "Tumor_Barcode", all.x = TRUE)
write.table(merged, file = "~/iCloud/dev/Sherlock-Lung/results2/log-reg-table4.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Optionally save results
#write.csv(results_df, "wGII_results.csv", row.names = FALSE)

#plot distribution of 


