library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)
library(EnhancedVolcano)
library(ggrepel)
library(ggvenn)
library(ggstatsplot)
library(ggplot2)
library(ggpubr)
library(sjPlot)
library(forcats)
library(ggsignif)
library(cowplot)
library(patchwork)
library(gridExtra)
library(grid)
library(VennDiagram)
library(hrbrthemes)
library(scales)
library(tidyr)
library(broom)
library(data.table)
library(flextable)
# Read the data
df <- read.csv("/Users/azhark/iCloud/dev/Sherlock-Lung/results2/deconvolution/Sherlock1+2+TCGA-LUAD+EAGLE_edge_normalized_logCPM_Danaher_mean.txt", sep = "\t", header = TRUE)

# Extract the first column to use as the index
index <- df[, 1]

# Remove the first column from the df
df2 <- df[, -1]

# Set row names from the extracted index
rownames(df2) <- index

# Optionally, transpose the data frame if needed
df2_transposed <- t(df2)
df3 <- data.frame(df2_transposed)

# Print the first few rows of the final data frame
head(df2_transposed)

setDT(df3, keep.rownames = "Sample")

#write out to tsv
write.table(df3, file = "~/iCloud/dev/Sherlock-Lung/results2/deconvolution/Sherlock1+2+TCGA-LUAD+EAGLE_edge_normalized_logCPM_Danaher_mean_transposed.txt", sep = "\t", quote = FALSE, row.names = FALSE)

deconv <- read.csv("~/iCloud/dev/Sherlock-Lung/results2/deconvolution/Sherlock1+2+TCGA-LUAD+EAGLE_edge_normalized_logCPM_Danaher_mean_transposed.txt", header = TRUE, sep="\t")

metadata = read.csv("~/iCloud/dev/Sherlock-Lung/data2/Sherlock_RNA-seq_batch1+2+EAGLE+TCGA-LUAD+TCGA-LUSC_all_samples_tmp.txt", sep = '\t', header=TRUE)
t = metadata %>% filter(Type == "Tumor")
t <- t[!duplicated(t$"RNAseq_SampleID"), ]
valid_rows <- t$RNAseq_SampleID[t$RNAseq_SampleID %in% deconv$Sample]

#subset rows for valid rows
deconv2 <- deconv[deconv$Sample %in% valid_rows,]

#read in ecDNA data
ecDNA <- read.csv("~/iCloud/dev/Sherlock-Lung/results2/ecDNA-annotations.txt", header = TRUE, sep="\t")
metadata2 = dplyr::select(metadata, WGS_Barcode, RNAseq_SampleID)
data <- merge(metadata2, ecDNA, by.x="WGS_Barcode", by.y="Tumor_Barcode")
merged2 <- merge(deconv2, data, by.x="Sample", by.y="RNAseq_SampleID")

# Transform the data into long format for easier plot
long_data <- merged2 %>%
  dplyr::select(ecDNA_status, Cell_Type, starts_with("x")) %>%
  pivot_longer(cols = -c(ecDNA_status, Cell_Type), names_to = "Sample", values_to = "Proportion") %>%
  mutate(Sample = gsub("x", "", Sample))  # Remove the "x" prefix from Sample names
#all the column sums add up to 1

custom_palette <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")

# Ensure that ecDNA_status and Cell_Type are factors
long_data$ecDNA_status <- as.factor(long_data$ecDNA_status)
long_data$Cell_Type <- as.factor(long_data$Cell_Type)

# Summarize data by averaging proportions for each ecDNA_status and Cell_Type
summarized_data <- long_data %>%
  group_by(ecDNA_status, Cell_Type) %>%
  summarize(Proportion = mean(Proportion))  # Use mean instead of sum

# Plot using the averaged proportions
ggplot(summarized_data, aes(x = ecDNA_status, y = Proportion, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "fill") +  # Use fill to normalize proportions within each ecDNA group
  scale_y_continuous(labels = scales::percent) +    # Convert y-axis to percentage
  scale_fill_manual(values = custom_palette) +      # Apply custom color palette
  labs(
    x = "ecDNA Status",
    y = "Proportion (as percentage)",
    fill = "Cell Type",
    title = "Proportion of Cell Types by ecDNA Status"
  ) +
  theme_minimal() +
  geom_text(aes(label = scales::percent(Proportion, accuracy = 1)),
            position = position_stack(vjust = 0.5), color = "white", size = 3)

# List of unique cell types
cell_types <- unique(long_data$Cell_Type)

# Create an empty list to store the plots
plot_list <- list()

# Loop over each cell type and create individual plots
for (cell in cell_types) {
  # Filter the data for the current cell type
  filtered_data <- long_data %>% filter(Cell_Type == cell)
  
  # Calculate the p-value using a Wilcoxon test (non-parametric test)
  wilcox_test <- wilcox.test(Proportion ~ ecDNA_status, data = filtered_data)
  p_value <- wilcox_test$p.value
  
  # Create the plot using ggbetweenstats
  plot <- ggbetweenstats(
    data = filtered_data,
    x = ecDNA_status,    # ecDNA status on the x-axis (comparing ecDNA+ vs ecDNA-)
    y = Proportion,      # Proportion on the y-axis
    type = "np",         # Non-parametric test (Wilcoxon or Kruskal-Wallis)
    pairwise.display = "significant",  # Only display significant pairwise comparisons
    pairwise.annotation = "p.value",   # Show p-values on the plot
    mean.plotting = FALSE,  # Show mean values on the plot
    median.plotting = FALSE, # Show mean values on the plot
    p.adjust.method = "holm",  # Adjust p-values using the Holm method
    palette = "Set2",          # Color palette
    results.subtitle = FALSE    # Display the main test result (e.g., Wilcoxon)
  ) + 
    labs(title = paste("Cell Type:", cell)) + # Add cell type as title
    theme_minimal() + 
    annotate("text", x = 1.5, y = max(filtered_data$Proportion) * 1.1, 
             label = paste0("p = ", signif(p_value, 3)), size = 5, color = "black")
  
    
  # Append the plot to the list
  plot_list[[cell]] <- plot
}

# Combine the individual plots using patchwork
combined_plot <- patchwork::wrap_plots(plot_list, ncol = 2)

# Display the combined plot
combined_plot

#columns are samples, rows are genes
met = read.delim('~/iCloud/dev/Sherlock-Lung/metadata/sherlock_1217_information_2023DEC21.tsv')
pge <- read.delim('~/iCloud/dev/Sherlock-Lung/results2/PGE.txt')
#count # of duplicates in Tumor_Barcode column of PGE
sum(duplicated(pge$Tumor_Barcode))



# Recodify histology
met$Histology = ifelse(met$Histology == 'Adenosquamous carcinoma',
                       'Other', met$Histology)
met$Histology = ifelse(met$Histology == 'Others',
                       'Other', met$Histology)

# Recodify Ancestry
met$Ancestry = ifelse(met$Assigned_Population %in% 
                        c('AFR', 'AMR or Mixed'),
                      'Other', met$Assigned_Population)

# Recodigy Stage
met$Stage_group = NA
met$Stage_group[met$Stage %in% c('I','IA','IA1',
                                 'IA2','IA3','IB')] = 'I'
met$Stage_group[met$Stage %in% c('II',  'IIA',  'IIB')] = 'II'
met$Stage_group[met$Stage %in% c('III', 'IIIA', 'IIIB')] = 'III'
met$Stage_group[met$Stage %in% c('IV',  'IVA')] = 'IV'

covdata = met %>% dplyr::select(Tumor_Barcode, Gender, Tumor_Purity,
                         Age, Histology, Stage_group, Ancestry, Smoking)

covdata$Gender = relevel(factor(covdata$Gender), ref = 'Female')
covdata$Histology = relevel(factor(covdata$Histology), ref = 'Adenocarcinoma')
covdata$Stage_group = relevel(factor(covdata$Stage_group), ref = 'I')
covdata$Ancestry = relevel(factor(covdata$Ancestry), ref = 'EUR')

merged3 <- merge(merged2, covdata, by.y="Tumor_Barcode", by.x="WGS_Barcode", all.x=TRUE)
merged4 <- merge(merged3, pge, by.y="Tumor_Barcode", by.x="WGS_Barcode", all.x=TRUE)

#replace NaN with 0 in PGE column 
merged4$ecDNA_bp[is.na(merged4$ecDNA_bp)] <- 0
#select all rows where PGE is > 0
merged4 <- merged4[merged4$ecDNA_bp > 0,]

#remove row with smoking unknown 
merged3 <- merged3[!is.na(merged3$Smoking),]

#subset for nonsmoker
merged3 <- merged3[merged3$Smoking == 'Non-Smoker',]
#subset for adeno
merged3 <- merged3[merged3$Histology == 'Adenocarcinoma',]

#run logistic regression on cell types
run_logistic_regression <- function(cell_type_column, data) {
    # Fit the logistic regression model
    model <- glm( ecDNA_status ~  data[[cell_type_column]] + Gender + Age + Ancestry + Tumor_Purity,# + Smoking, 
                 family = binomial, data = data)
    
    # Check for perfect separation
    if (summary(model)$deviance == 0) {
      warning(paste("Model fitting failed for", cell_type, "due to perfect separation"))
      return(NULL)
    }
    
    # Extract model summary
    tidy_model <- tidy(model)
    
    # Calculate confidence intervals and odds ratios
    conf_intervals <- exp(confint(model))
    tidy_model <- tidy_model %>%
      mutate(conf.low = conf_intervals[, 1],
             conf.high = conf_intervals[, 2],
             odds_ratio = exp(estimate))
    
    # Add gene name to the results
    tidy_model$cell_type <- cell_type_column
    
    return(tidy_model)
}

# Initialize results list
results_list <- list()

merged3 = merged3 %>% filter(Smoking != "Unknown")
smokers = merged3 %>% filter(Smoking == 'Smoker')
nonsmokers = merged3 %>% filter(Smoking == 'Non-Smoker')


cell_types <- c("B.cells", "CD45", "CD8.T.cells", "Cytotoxic.cells", 
                                "DC", "Exhausted.CD8", "Macrophages", "Mast.cells", 
                                "Neutrophils", "NK.CD56dim.cells", "NK.cells", "T.cells", 
                                "Th1.cells", "Treg")
for (cell in cell_types) {
  print(paste("Running logistic regression for gene:", cell))
  model_summary <- run_logistic_regression(cell, merged3)
  if (!is.null(model_summary)) {
    results_list[[cell]] <- model_summary
    print(paste("Completed logistic regression for gene:", cell))
  }
}
  
# Combine all results into a single dataframe
results_summary <- bind_rows(results_list)

#remove all the rows with the term intercept
df_cleaned <- results_summary %>% filter(term == "data[[cell_type_column]]")

# Adjust p-values using FDR
df_cleaned$fdr <- p.adjust(df_cleaned$p.value, method = "fdr")

# Sort by lowest FDR p-value first, and then highest odds ratio
df_cleaned <- df_cleaned %>%
  arrange(fdr, desc(odds_ratio))

# The final results_summary dataframe contains the summarized logistic regression results
print(df_cleaned)

df_cleaned$term <- df_cleaned$cell_type

df_cleaned %>% 
  ggplot(aes(log2(odds_ratio),-log10(fdr),fill=cell_type))+
  geom_hline(yintercept = -log10(0.05),col="red",linetype=2)+
  geom_hline(yintercept = -log10(0.01),col="red",linetype=2)+
  geom_vline(xintercept = 0,col="gray10",linetype=1,linewidth=0.5)+
  geom_point(pch=21,size=3.5,col='black',stroke=0.2)+
  ggrepel::geom_text_repel(data=df_cleaned %>% filter(fdr<1),aes(label=cell_type,col=cell_type),force=20)+
  scale_x_continuous(breaks = pretty_breaks(n = 7),limits = c(-0.5,0.5))+
  scale_y_continuous(breaks = pretty_breaks(n=7))+
  theme_bw(base_size = 12)+ # Changed from theme_ipsum_rc to theme_bw
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5)
  )+
  labs(x = 'Odd ratio (log2)', y = '-log10(FDR)',title='Immune cell composition')+
  guides(fill="none",color='none')+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  coord_cartesian(clip = 'off')

# Round all columns except for fdr to 2 decimal places
df_cleaned <- df_cleaned %>%
  mutate(
    across(c(estimate, std.error, statistic, p.value, conf.low, conf.high, odds_ratio), round, 2),
    fdr = formatC(fdr, format = "e", digits = 2)
  )



# Assuming your data is in a data frame called df
# Create a publication-quality table
flex_table <- flextable(df_cleaned)$ecDNA_bp
flex_table <- set_caption(flex_table, "Association between immune cell types and ecDNA, corrected for age. gender, histology, ancestry, and tumor purity")

flex_table <- flex_table %>%
  set_table_properties(width = .8, layout = "autofit") %>%
  theme_zebra() %>%
  align_text(align = "center", part = "all") %>%
  align(align = "center", part = "header") %>%
  autofit()

# Save the table as a Word document
save_as_docx(flex_table, path = "/Users/khandekara2/Library/Mobile Documents/com~apple~CloudDocs/Documents/Docs/ecDNA_manuscript/Supplementary/immune_supp_table.docx")

# Save the table as a PowerPoint slide
save_as_pptx(flex_table, path = "/Users/khandekara2/Library/Mobile Documents/com~apple~CloudDocs/Documents/Docs/ecDNA_manuscript/Supplementary/immune_supp_figure.pptx")


