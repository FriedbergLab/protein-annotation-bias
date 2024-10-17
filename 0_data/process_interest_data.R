library(tidyverse)

# Helper functions
percentage_increase <- function(row) {
  numeric_cols <- row[sapply(row, is.numeric)]
  first_non_zero <- numeric_cols[numeric_cols != 0][1]
  if (!is.numeric(first_non_zero)) {
    return(NA) # If there is no non-zero and numeric value, return NA
  }
  last_value <- as.numeric(row[length(numeric_cols)])
  # Calculate the percentage increase
  percentage_increase <- round((last_value / first_non_zero),4)
  return(percentage_increase)
}

absolute_increase <- function(row) {
  numeric_cols <- row[sapply(row, is.numeric)]
  first_non_zero <- numeric_cols[numeric_cols != 0][1]
  if (!is.numeric(first_non_zero)) {
    return(NA) # If there is no non-zero and numeric value, return NA
  }
  last_value <- as.numeric(row[length(numeric_cols)])
  # Calculate the percentage increase
  absolute_increase <- round((last_value - first_non_zero),4)
  return(absolute_increase)
}

# GO_based data
proteins_with_IC <- read.csv(file="0_data/processed_GO/AllProteins_2022.tsv",sep="\t",header=T)%>%
  select(Genes,Count)%>%
  filter(Count > 0)%>%
  select(Genes)%>%
  unlist()%>%
  unique()

# TODO: maybe only use STRING, or add Biopython parser script
# ID mapping from UniProt-SwissProt: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/, parsed using Biopython # 
proteinstring <- read.csv(file="0_data/ext_data/sw-human-2024.txt",sep="\t",header=F)%>%
  select(2,3,5)%>%
  rename_with(~c("Genes","Common","STRING"))
# Additional STRING aliases from STRING: https://stringdb-downloads.org/download/protein.aliases.v12.0.txt.gz
proteinstring_add <- read.csv(file="0_data/ext_data/9606.protein.aliases.v12.0.txt",sep="\t",header=T)%>%
  rename_with(~c("STRING","alias","source"))%>%
  mutate(STRING = gsub(".*\\.(ENSP[0-9]+)", "\\1",STRING))%>%
  filter(STRING != alias)

## Number of FPEs from fractional counts ####
# Around 18,000 proteins
proteincounts <- read.csv(file="0_data/ext_data/protein_counts.tsv",sep="\t",header=F)%>%
  rename_with(~c("STRING","Date","FractionalCount"))%>%
  filter(Date < 2023)%>% # cumulative sum everything before 2023
  group_by(STRING)%>%
  arrange((Date))%>%
  mutate(FractionalCount2 = cumsum(FractionalCount))%>%select(-3)%>%ungroup()%>%
  complete(STRING,Date,fill=list(FractionalCount2=NA))%>%
  group_by(STRING)%>%
  fill(FractionalCount2,.direction = "down")%>%
  filter(Date >= 2013)
fc_data <- proteincounts%>%
  pivot_wider(names_from = "Date",values_from="FractionalCount2")%>%
  left_join(proteinstring,by="STRING")%>%
  select(1,12,13,2:11)%>%
  ungroup()%>%
  mutate_at(vars(4:13), ~replace(., is.na(.), 0))

# For proteins that do have counts but their STRING is different; only match with proteins that have IC
alias_withIC <- proteinstring_add%>%
  filter((alias %in% proteins_with_IC) & alias %in% proteinstring$Genes)%>%
  distinct(STRING, alias, .keep_all = T)%>%
  rename(Genes = alias)%>%
  select(-source)
missing_fc_data <- fc_data%>%
  filter(is.na(Genes) & `2022` > 0)%>%
  select(-c(Genes,Common))%>%
  left_join(alias_withIC,by="STRING")
notmissing_fc_data <- missing_fc_data%>%
  filter(!is.na(Genes))%>%
  select(1,12,2:11)%>%
  rbind(fc_data%>%select(-Common)%>%filter(!is.na(Genes)))%>%
  group_by(Genes)%>%
  filter(`2022` == max(`2022`))%>%
  ungroup()
stillmissing <- missing_fc_data%>%
  filter(is.na(Genes))%>%select(-Genes)%>%
  rename(alias = STRING)%>%
  left_join(proteinstring_add,by="alias")%>%
  select(-source)%>%
  left_join(alias_withIC,by="STRING")%>%
  filter(!is.na(Genes))%>% # these proteins do not have known annotations in selected evidence 
  distinct(.)%>%
  select(12,13,2:11)
filtered_fc <- rbind(stillmissing, notmissing_fc_data)
filtered_fc <- filtered_fc%>%
  group_by(Genes)%>%
  filter(`2022` == max(`2022`))%>%
  ungroup()%>%group_by(STRING)%>%
  filter(`2022` == max(`2022`))%>%
  ungroup()%>%
  distinct(STRING, Genes, .keep_all = T) # Some STRING aliases and UniProt duplicates are still in, filter using DISEASES and IC later

# Removing clinical biomarkers ####
# Diagnostic protein marker names from MarkerDB: https://markerdb.ca/pages/download_all_diagnostic_proteins?format=tsv
markerDB_names <- read.csv(file="0_data/ext_data/all_diagnostic_proteins.tsv",header=T,sep="\t")%>%
  select(c("name","gene_name","uniprot_id"))
# Diagnostic protein markers from MarkerDB: https://markerdb.ca/pages/download_all_proteins?format=tsv
markerDB <- read.csv(file="0_data/ext_data/all_proteins.tsv",header=F,sep="\t")%>%
  select(V2)%>%
  rename(name = V2)%>%
  full_join(markerDB_names,by="name")%>%
  filter(uniprot_id != "" | gene_name != "")%>%
  distinct(.)%>%
  rename(Common = gene_name)%>%
  left_join(proteinstring,by="Common")%>%
  mutate(uniprot_id = ifelse(uniprot_id == "",Genes, uniprot_id))%>%
  mutate(Genes = ifelse(is.na(Genes),uniprot_id,Genes))
# remove_markerDB <- markerDB%>%
#   filter(Genes %in% filtered_fc$Genes)%>%
#   left_join(filtered_fc%>%select(-STRING),by="Genes")

# With clinical proteins
clinical_fc <- filtered_fc%>%
  left_join(proteinstring%>%select(-STRING),by="Genes")%>%
  select(1,2,13,3:12)
# Without clinical proteins
fc <- filtered_fc%>%
  filter(!(Genes %in% markerDB$Genes))%>%
  filter(!(Genes %in% c("P01730", "P01375")))%>% # remove TNF and CD4, not included in MarkerDB
  left_join(proteinstring%>%select(-STRING),by="Genes")%>%
  select(1,2,13,3:12)
fc = fc%>%
  ungroup()%>%
  rowwise()%>%mutate(rel_gain2=percentage_increase(c_across(4:13)))%>% # percentage increase is 2022 value divided by first non-zero value since 2013
  rename_with(~c(colnames(fc_data),"rel_gain"))%>%
  rowwise()%>%mutate(gain=absolute_increase(c_across(4:13)))%>%
  mutate(gain = ifelse((is.na(gain) & `2013` == 0 & `2022` == 0), 0, gain), rel_gain = ifelse((is.na(rel_gain) & `2013` == 0 & `2022` == 0), 0, rel_gain))
# Duplicated UniProt entries with different STRING
fc%>%group_by(Genes)%>%summarise(n=n())%>%filter(n>1L)
# Duplicated STRING entries with different UniProt
fc%>%group_by(STRING)%>%summarise(n=n())%>%filter(n>1L)

saveRDS(clinical_fc, file="1_metrics/interest_FPE_clinical.RDS")
saveRDS(fc, file="1_metrics/interest_FPE.RDS")
