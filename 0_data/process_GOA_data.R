library(tidyverse)

process_protein_files <- function(folder_path, identifier) {
  # Get list of files matching the pattern
  files <- list.files(folder_path, pattern = paste0(identifier,"Proteins_\\d{4}\\.tsv"), full.names = TRUE)
  # Read all files and combine them into a single dataframe
  all_data <- map_df(files, ~{
    year <- str_extract(.x, "\\d{4}")
    read_tsv(.x) %>% mutate(Date = year)
  })
  # Separate dataframes by Aspect
  aspects <- unique(all_data$Aspect)
  
  # Process each aspect
  result <- map(aspects, function(asp) {
    all_data %>%
      filter(Aspect == asp) %>%
      select(-Aspect) %>%
      pivot_longer(cols = c(Count, UniCount, IA),
                   names_to = "Metric",
                   values_to = "Value") %>%
      pivot_wider(names_from = Date,
                  values_from = Value,
                  values_fill = list(Value = 0)) %>%
      pivot_longer(cols = -c(Genes, Metric),
                   names_to = "Date",
                   values_to = "Value") %>%
      pivot_wider(names_from = Metric,
                  values_from = Value)
  }) %>% setNames(aspects)
  
  return(result)
}

for (identifier in c("All", "HTP", "LTP")) {
  combined_data <- process_protein_files(folder_path = "0_data/processed_GO", identifier = identifier)
  saveRDS(combined_data, paste0("1_metrics/GO_exp_aspect_", identifier, ".RDS"))
}
