# Genomic epidemiology of hospital related infections in Hong Kong 2022

This project focuses on analyzing hospital-related infections in Hong Kong, including case classifications, genomic surveillance, and mobility data. The repository contains scripts, data, and results for visualizing and interpreting the findings.

## Project Structure

### 1. **Data**
The `data/` directory contains raw and processed data files used in the analysis. Key subdirectories include:
- **`case_curve/`**: Contains JSON and Excel files for aggregated case data.
  - `HK_case_data.json`: JSON file with aggregated confirmed case classifications.
  - `hk_epicurve.xlsx`: Processed case data in Excel format.
- **`hospital_data/`**: Contains metadata for hospitals and districts.
  - `metadata_cases_plot.csv`: Metadata for hospital-acquired infections.
- **`mobility_data/`**: Contains mobility and transport data.
  - `2022_HK_Region_Mobility_Report.csv`: Google Community Mobility Report for Hong Kong.
  - `statistics_on_daily_passenger_traffic.csv`: Cross-border passenger traffic data.
  - `public_transport.xlsx`: Local transport data.

### 2. **Scripts**
The `scripts/` directory contains R scripts for data processing and visualization. Key scripts include:
1. **metadata_visualization.R**  
   - Reads and processes metadata from various sources.  
   - Generates figures and summarizing tables about case data, such as Panel A/B epidemiological curves.  

2. **tree_visualization.R**  
   - Loads and processes MCC (BEAST) and ML phylogenetic trees.  
   - Annotates hospital clusters, collapses large community clades, and creates publication-quality tree figures.  

3. **GAM.R**  
   - Reads Markov jump data and predictor datasets.  
   - Builds and diagnoses a GAM (Negative Binomial) model, performs VIF checks, and creates time-series plots of residual autocorrelation.  


### 3. **Results**
The `results/` directory contains output figures and processed data

## How to Run
1. **Install Required Libraries**:
   Ensure the following R libraries are installed:
   - `tidyverse`, `jsonlite`, `lubridate`, `ggplot2`, `ggforce`, `ggrepel`, `MetBrewer`, `ggbump`, `scales`, `sp`, `gridExtra`, `RColorBrewer`, `readxl`, `writexl`.

2. **Run Scripts**:
   - Execute `1.metadata_visualization.R` to generate Figures 1a, 1b, 1c, and 1d.
   - Execute `2.tree_visualization.R` to generate phylogenetic tree visualizations.
   - Execute `3.GAM.R` to generate GAM results.
   - Execute `GISAID_data_processing/QC.R` to process and filter GISAID genomic data. (Details refer to [here](scripts/GISAID_data_processing/README.md))

3. **View Results**:
   - Figures are saved in the `results/` directory.

## Citation
Pending

---