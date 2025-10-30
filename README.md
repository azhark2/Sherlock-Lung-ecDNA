This repository contains the code required to reproduce all the results of 'Examining the Role of Extrachromosomal DNA in 1,216 Lung Cancers'
(https://www.biorxiv.org/content/10.1101/2025.06.03.657117v2).

## 🚀 Quick Start

### Prerequisites
- R version 4.0 or higher
- Required R packages (see Installation section)

### Installation
```r
# Install required packages
required_packages <- c(
  "sjPlot", "sjlabelled", "sjmisc", "ggplot2", "EnhancedVolcano",
  "ggrepel", "ggvenn", "ggstatsplot", "ggpubr", "forcats",
  "ggsignif", "cowplot", "patchwork", "gridExtra", "grid",
  "VennDiagram", "hrbrthemes", "scales", "tidyr", "broom",
  "data.table", "flextable", "tidyverse", "stringr", "dplyr",
  "ggalluvial", "logistf", "openxlsx", "survival", "survminer",
  "nlstools", "UpSetR", "rstatix"
)
```


### Install missing packages
```r
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
```

### Analysis Overview
1. **Main_Figures.Rmd:** Main manuscript figures
2. **Supplementary_Figures.Rmd:** Extended analyses and supplementary visualizations
3. **Supplementary_Tables.Rmd:** Statistical tables and data


**Data files are not included in this repository.** Set up this directory structure:

```
data/
├── metadata/
│   ├── sherlock_1217_information_2023DEC21.tsv
│   ├── Country_group_for_Azhar_2024-12-13.csv
│   └── [other metadata files]
└── results/
    ├── log-reg-table5.tsv
    ├── log-reg-results2/
    └── [other result files]
```

### Data Availability
Contact the first author for data access (azhark@gmail.com)

### Running the Analysis
1. Set up the data directory structure as described above
2. Open R/RStudio and set working directory to the repository root
3. Run analysis scripts from the scripts/ directory


