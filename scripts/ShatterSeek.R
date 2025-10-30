#script to run ShatterSeek on Sherlock data
library(ShatterSeek)
library(gridExtra)
library(tidyverse)
library(scales)
library(ggrepel)
library(broom)
library(GenomicRanges)
library(gridExtra)
library(cowplot)
data(DO17373)

# CN_data <- CNVsegs(chrom=as.character(SCNA_DO17373$chromosome),
#                    start=SCNA_DO17373$start,
#                    end=SCNA_DO17373$end,
#                    total_cn=SCNA_DO17373$total_cn)
# 
# SV_data <- SVs(chrom1=as.character(SV_DO17373$chrom1),
#                pos1=as.numeric(SV_DO17373$start1),
#                chrom2=as.character(SV_DO17373$chrom2),
#                pos2=as.numeric(SV_DO17373$end2),
#                SVtype=as.character(SV_DO17373$svclass),
#                strand1=as.character(SV_DO17373$strand1),
#                strand2=as.character(SV_DO17373$strand2))
# 
# # Run ShatterSeek
# chromothripsis <- shatterseek(
#   SV.sample=SV_data,
#   seg.sample=CN_data,
#   genome="hg38")

#you'll need to use the CNVkit .cns files (Biowulf: /data/khandekara2/dev/Sherlock-Lung/data/cnvkit_corrected_seeds), not the Battenberg ones 

#change directory
setwd("~/iCloud/dev/Sherlock-Lung/data/ShatterSeek-input/")
CN_files <- list.files(path = "~/iCloud/dev/Sherlock-Lung/data/ShatterSeek-input/CNVkit/", pattern = ".segments.tsv", full.names = FALSE)
SV_files <- list.files(path = "~/iCloud/dev/Sherlock-Lung/data/ShatterSeek-input/", pattern = ".SV.bedpe.tsv", full.names = FALSE)

#sort lists aplhabetically
CN_files <- sort(CN_files)
SV_files <- sort(SV_files)

# Initialize an empty list to store the results for each sample
results_list <- list()

# Initialize counters for criteria not met
criteria_counters <- list(
  pval_fragment_joins_not_met = 0,
  chr_breakpoint_enrichment_not_met = 0,
  pval_exp_cluster_not_met = 0,
  number_events_not_met = 0,
  number_TRA_not_met = 0,
  oscillating_CN_segments_not_met = 0
)

all_chrom_summaries <- list()


# Create an empty list to store plots
plot_list <- list()
bb_dir <- "/Users/azhark/iCloud/dev/Sherlock-Lung/data/CN/bed/"
#for loop to iterate through all samples in directory and run ShatterSeek
for (i in 1:length(SV_files)) {
  #read in structural variant data
  SV <- read.table(SV_files[i], sep = "\t", header = TRUE)
  
  sample_name <- gsub(".SV.bedpe.tsv", "", basename(SV_files[i]))
  print(sample_name)
  print(paste("Processing:", paste0(sample_name, ".segments.tsv"), "and", SV_files[i]))
  
  CN_file <- paste0(sample_name, ".segments.tsv")
  
  #read in copy number data if it exists
  # Check if the CN file exists
  if (!file.exists(CN_file)) {
    print(paste("File not found:", CN_file, "- Skipping to next sample."))
    next  # Skip to the next iteration of the loop if the file does not exist
  }
  
  CN <- read.table(CN_file, sep = "\t", header = TRUE)

  #MERGE SEGMENTS
  # d is a data.frame with colums: chr, start, end, total_cn
  dd <- CN
  dd$total_cn[dd$total_cn == 0] <- 150000
  dd$total_cn[is.na(dd$total_cn)] <- 0
  dd <- as(dd,"GRanges")
  cov <- coverage(dd,weight = dd$total_cn)
  dd1 <- as(cov,"GRanges")
  dd1 <- as.data.frame(dd1)
  dd1 <- dd1[dd1$score !=0,]
  dd1 = dd1[,c(1,2,3,6)]
  names(dd1) <- names(CN)[1:4]
  dd1$total_cn[dd1$total_cn == 150000] <- 0
  CN= dd1; rm(dd)
  

  # Ensure that chromosome names are in the expected format (numeric or "X")
  valid_chromosomes <- as.character(c(1:22, "X"))
  CN <- CN[CN$chr %in% valid_chromosomes, ]
  SV <- SV[SV$chrom1 %in% valid_chromosomes & SV$chrom2 %in% valid_chromosomes, ]

  # Get sample name from the CN file name
  sample_name <- strsplit(basename(CN_files[i]), ".segments.tsv")[[1]][1]
  
  #run ShatterSeek
  CN_data <- CNVsegs(chrom=as.character(CN$chr),
                     start=CN$start,
                     end=CN$end,
                     total_cn=CN$total_cn)
  
  SV_data <- SVs(chrom1=as.character(SV$chrom1), 
                 pos1=as.numeric(SV$start1),
                 chrom2=as.character(SV$chrom2), 
                 pos2=as.numeric(SV$end2),
                 SVtype=as.character(SV$svclass), 
                 strand1=as.character(SV$strand1),
                 strand2=as.character(SV$strand2))
  
  ecDNAdir <- "/Users/azhark/iCloud/dev/Sherlock-Lung/results2/ecDNA-regions/rescaled/temp/"
  ecDNA_detected <- FALSE
  if (file.exists(paste0(ecDNAdir, sample_name, ".bed"))) {
    ecDNA_detected <- TRUE
    ecDNA_file <- read.table(paste0(ecDNAdir, sample_name, ".bed"), sep = "\t", header = FALSE)
    #get all the unique values in the first column, an strip the 'chr' prefix if it exists
    ecDNA_chr <- unique(ecDNA_file[, c(1)])
    ecDNA_chr <- gsub("chr", "", ecDNA_chr)
    print(ecDNA_chr)
  }
  
  chromothripsis <- shatterseek(
    SV.sample=SV_data,
    seg.sample=CN_data,
    genome="hg38")
  
  chromSummary_df <- chromothripsis@chromSummary
  chromSummary_df$Sample <- sample_name
  all_chrom_summaries[[sample_name]] <- chromSummary_df
  #save chromSummary_df
  #write.table(chromSummary_df, file = paste0("/Users/khandekara2/iCloud/dev/Sherlock-Lung/data/ShatterSeek-input/output/", sample_name, "_chromSummary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Initialize a vector to store chromosomes with chromothripsis
  chrom_with_chromothripsis <- c()
  
  #High confidence: at least 6 interleaved intrachromosomal SVs, 7 contiguous segments oscillating between 2 CN states, 
  #the fragment joins test, and either the chromosomal enrichment or the exponential distribution of breakpoints test.
  
  #High confidence: at least 3 interleaved intrachromosomal SVs and 4 or more interchromosomal SVs, 7 contiguous segments oscillating between 2 CN states and the fragment joins test.
  
  #Low confidence: at least 6 interleaved intrachromosomal SVs, 4, 5 or 6 adjacent segments oscillating between 2 CN states, the fragment joins test, and either the chromosomal enrichment or the exponential distribution of breakpoints test.
  
  # Loop through each row in chromSummary_df (each chromosome)
  for (j in 1:nrow(chromSummary_df)) {
    # Check each criterion and update the counters if not met
    pval_fragment_joins_met <- !is.na(chromSummary_df$pval_fragment_joins[j]) && chromSummary_df$pval_fragment_joins[j] > 0.05
    if (!pval_fragment_joins_met) criteria_counters$pval_fragment_joins_not_met <- criteria_counters$pval_fragment_joins_not_met + 1
    
    # Modify the criteria to satisfy if either chr_breakpoint_enrichment < 0.05 OR pval_exp_cluster < 0.05
    chr_breakpoint_or_exp_cluster_met <- (!is.na(chromSummary_df$chr_breakpoint_enrichment[j]) && chromSummary_df$chr_breakpoint_enrichment[j] < 0.05) ||
      (!is.na(chromSummary_df$pval_exp_cluster[j]) && chromSummary_df$pval_exp_cluster[j] < 0.05)
    
    if (!chr_breakpoint_or_exp_cluster_met) criteria_counters$chr_breakpoint_enrichment_not_met <- criteria_counters$chr_breakpoint_enrichment_not_met + 1
    
    number_events_met <- !is.na(chromSummary_df$number_DEL[j]) && !is.na(chromSummary_df$number_DUP[j]) &&
      !is.na(chromSummary_df$number_h2hINV[j]) && !is.na(chromSummary_df$number_t2tINV[j]) &&
      (chromSummary_df$number_DEL[j] + chromSummary_df$number_DUP[j] + 
         chromSummary_df$number_h2hINV[j] + chromSummary_df$number_t2tINV[j]) >= 6
    if (!number_events_met) criteria_counters$number_events_not_met <- criteria_counters$number_events_not_met + 1
    
    number_TRA_met <- !is.na(chromSummary_df$number_TRA[j]) && chromSummary_df$number_TRA[j] >= 4
    if (!number_TRA_met) criteria_counters$number_TRA_not_met <- criteria_counters$number_TRA_not_met + 1
    
    oscillating_CN_segments_met <- !is.na(chromSummary_df$max_number_oscillating_CN_segments_2_states[j]) &&
      chromSummary_df$max_number_oscillating_CN_segments_2_states[j] >= 7
    if (!oscillating_CN_segments_met) criteria_counters$oscillating_CN_segments_not_met <- criteria_counters$oscillating_CN_segments_not_met + 1
    
    # Determine if all criteria are met for this chromosome
    criteria_met <- pval_fragment_joins_met && chr_breakpoint_or_exp_cluster_met && number_events_met && number_TRA_met && oscillating_CN_segments_met
    #print(criteria_met, pval_fragment_joins_met, chr_breakpoint_or_exp_cluster_met, number_events_met, number_TRA_met, oscillating_CN_segments_met)
    
    # If criteria are met for this chromosome, add it to the list
    if (!is.na(criteria_met) && criteria_met) {
      chrom_with_chromothripsis <- c(chrom_with_chromothripsis, as.character(chromSummary_df$chromosome[j]))
    }
    
    #PLOTTING
    if (ecDNA_detected && !criteria_met && !is.null(ecDNA_chr)) {
      chr_value <- as.character(chromSummary_df$chrom[j])
      if (chr_value %in% ecDNA_chr) {
        print(paste("Attempting to create plot for:", sample_name, "chromosome:", chr_value))
        
        chromothripsis_plot <- tryCatch({
          plots <- plot_chromothripsis(
            ShatterSeek_output = chromothripsis,
            chr = chr_value,
            sample_name = sample_name,
            genome = "hg38"
          )
          print("Plot list created successfully")
          plots  # This is a list of plot objects
        }, error = function(e) {
          print(paste("Error creating plot:", e$message))
          return(NULL)
        })
        
        # Store the plot list if created successfully
        if (!is.null(chromothripsis_plot)) {
          key <- paste0(sample_name, "_", chr_value)
          print(paste("Storing plot with key:", key))
          plot_list[[key]] <- chromothripsis_plot
          print(paste("Current number of plots in list:", length(plot_list)))
        }
      }
    }
    
    
  }
  
  chrom_with_chromothripsis <- chromSummary_df %>%
    dplyr::filter(
      pval_fragment_joins >= 0.05,
      (chr_breakpoint_enrichment < 0.05 | pval_exp_cluster < 0.05),
      (number_DEL + number_DUP + number_h2hINV + number_t2tINV) >= 3,
      number_TRA >= 4,
      max_number_oscillating_CN_segments_2_states >= 4
    ) %>%
    pull(chrom)
  
  chromothripsis_df <- chromSummary_df %>%
    dplyr::filter(
      pval_fragment_joins >= 0.05,
      (chr_breakpoint_enrichment < 0.05 | pval_exp_cluster < 0.05),
      (number_DEL + number_DUP + number_h2hINV + number_t2tINV) >= 3,
      number_TRA >= 4,
      max_number_oscillating_CN_segments_2_states >= 4
    )
  
  #remove entires that have NA in start or end columns
  chromothripsis_df <- chromothripsis_df[!is.na(chromothripsis_df$start) | !is.na(chromothripsis_df$end),]
  #now make a sample column that has the sample name, make it the same length as the chromothripsis_df
  chromothripsis_df$Sample <- rep(sample_name, nrow(chromothripsis_df))
  
  #replace all empty strings with NA
  chromothripsis_df <- chromothripsis_df %>% 
    mutate(across(everything(), ~ if_else(. == "", NA, .)))
  
  #subset for only chrom, start, end, and sample columns
  chromothripsis_df <- chromothripsis_df %>% 
    dplyr::select(chrom, start, end, Sample)
  
  #save chrom_with_chromothripsis
  print(paste0("/Users/khandekara2/iCloud/dev/Sherlock-Lung/results2/chromothripsis-regions/", sample_name, "_chromothripsis.tsv"))
  #write.table(chromothripsis_df, file = paste0("/Users/khandekara2/iCloud/dev/Sherlock-Lung/results2/chromothripsis-regions/", sample_name, "_chromothripsis.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  
  # Determine if chromothripsis was detected and concatenate the list of chromosomes
  chromothripsis_detected <- length(chrom_with_chromothripsis) > 0
  chrom_with_chromothripsis_str <- ifelse(chromothripsis_detected, paste(chrom_with_chromothripsis, collapse = ", "), "None")
  
  # Create a new data frame with the sample name, YES/NO, and the chromosomes with chromothripsis
  results_table <- data.frame(
    Sample = sample_name,  # Assign the sample name
    Chromothripsis = ifelse(chromothripsis_detected, "YES", "NO"),
    Chromosomes = chrom_with_chromothripsis_str
  )
  
  # Add the results_table to the list
  results_list[[i]] <- results_table
  

}

#print the statistical criteria that was used 

# After the main loop, combine all the dataframes:
final_chrom_summary <- do.call(rbind, all_chrom_summaries)

# If you want to reorder the columns to have Sample as the first column:
final_chrom_summary <- final_chrom_summary[, c("Sample", setdiff(names(final_chrom_summary), "Sample"))]

# Combine all data frames in the list into a single data frame
final_results <- do.call(rbind, results_list)

# View the final results
#print(final_results)

# Save the combined data frame to a CSV file
# write.csv(final_results, file = "combined_high_confidence_results.csv", row.names = FALSE)

# Summary statistics
summary_stats <- final_chrom_summary %>%
  group_by(Sample) %>%
  summarise(
    mean_events = mean(number_DEL + number_DUP + number_h2hINV + number_t2tINV, na.rm = TRUE),
    max_events = max(number_DEL + number_DUP + number_h2hINV + number_t2tINV, na.rm = TRUE),
    mean_oscillating_CN = mean(max_number_oscillating_CN_segments_2_states, na.rm = TRUE),
    min_pval_fragment_joins = min(pval_fragment_joins, na.rm = TRUE),
    min_chr_breakpoint_enrichment = min(chr_breakpoint_enrichment, na.rm = TRUE)
  )

print(table(final_results$Chromothripsis))

print(summary_stats)

#write out final results to tsv file
#write.table(final_results, file = "combined_high_confidence_results.tsv", sep = "\t", row.names = FALSE)

print(paste("Total number of plots to save:", length(plot_list)))
print("Plot keys:")
print(names(plot_list))

# Save plots to PDF with proper arrangement
pdf("/Users/azhark/iCloud/dev/Sherlock-Lung/figures/ecDNA-no-chromothripsis.pdf", width=12, height=15)

if (length(plot_list) > 0) {
  for (key in names(plot_list)) {
    print(paste("Saving plot:", key))
    plot_components <- plot_list[[key]]
    
    # Create a new page for each set of plots
    grid.newpage()
    
    # Arrange all components vertically
    # Calculate heights based on the number of components
    heights <- c(0.2, 0.3, 0.2, 0.3)  # Adjust these ratios as needed
    if (length(plot_components) == 5) {  # If BAF plot is included
      heights <- c(0.15, 0.25, 0.2, 0.2, 0.2)
    }
    
    # Create layout
    grid.arrange(
      grobs = plot_components,
      heights = heights,
      ncol = 1
    )
  }
} else {
  print("No plots found in plot_list!")
}

dev.off()

print("PDF creation completed")





