# COVID-19 Data Analysis in Hong Kong (2022)

The analysis conducted on the Hong Kong COVID-19 data covers the period from **January 1, 2022, to December 31, 2022**. The dataset includes **12,641 cases**.

## Data Filtering Process

1. **Quality Control**:
   - A total of **561 sequences** with poor quality were filtered out using **Nextclade** and **Pangolin**.

2. **Exclusion of Specific Samples**:
   - To exclude samples sequenced by **PolyU**, an additional **134 cases** were removed.

### Final Dataset

After the filtering process, the final dataset comprises **11,946 cases** of high-quality sequences.

---

## Steps to Process Data

```bash
Step 1: Run Quality Control Script
# Run the QC script
Rscript QC.R

Step 2: Extract FASTA Sequences
# Extract high-quality sequences using Python script
python ./Scripts/extract_fasta.py https://hkuhk-my.sharepoint.com/:u:/g/personal/limt_hku_hk/EbWRGwXzKTZOgh8HUW2KXkUBmnDbIo4ygRQkpQqHoeiWbA?e=6AvCoG ./Data/filter.list extract.fasta
