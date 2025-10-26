#########################step1: pangolin########################################
setwd("~/pangolin")
system("conda run -n pangolin pangolin --update")
file_current_fasta<-"https://hkuhk-my.sharepoint.com/:u:/r/personal/limt_hku_hk/Documents/Update_GISAID/gisaid.fasta?csf=1&web=1&e=B1HRZA"
system(paste0("conda run -n pangolin pangolin ", file_current_fasta))

#########################step2: nextclade#######################################
setwd("~/nextclade")
file_current_fasta <-"https://hkuhk-my.sharepoint.com/:u:/r/personal/limt_hku_hk/Documents/Update_GISAID/gisaid.fasta?csf=1&web=1&e=B1HRZA"
system(paste0("./Scripts/run_nextclade.sh ", file_current_fasta, " ./"))

#########################step3: filter##########################################
###1.filter bad quality samples
library(writexl)
library(tidyverse)
library(Biostrings)
library(dplyr)
library(stringr)

df_lin_nextclade <- read_tsv("~/nextclade/nextclade.tsv")
df_lin_pangolearn <- read_csv("~/pangolin/lineage_report.csv")

id <- union(
  df_lin_nextclade[df_lin_nextclade$qc.overallStatus %in% c("bad", "mediocre"), ] %>% .$seqName,
  df_lin_pangolearn[df_lin_pangolearn$qc_status %in% c("fail"), ] %>% .$taxon
)

meta<-readxl::read_xlsx("~/Data/metadata_gisaid.xlsx")
meta <- meta %>% 
  `rownames<-`(.$strain) %>% 
  .[setdiff(.$strain, id), , drop = FALSE]

###2.filter PolyU samples
meta.polyu<-readxl::read_xlsx("~/Data/metadata_PolyU.xlsx")
meta.polyu<-tidyr::separate(meta.polyu,FASTA,into=c("id","other"),sep="/")

matches <- grepl(paste(meta.polyu$id, collapse="|"), meta$strain)
meta_filtered <- meta[!matches, ]

write_xlsx(meta_filtered,"~/Data/metadata_filtered.xlsx")






