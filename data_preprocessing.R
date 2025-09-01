# Load all necessary libraries
library(httr)
library(readr) 
library(dplyr)
library(tidyverse)
library(data.table)
library(jsonlite)
library(stringr)
library(purrr)
library(furrr)
library(UniProt.ws)
library(UniprotR)
library(Biostrings)
library(ape)
library(R.utils)
library(pwalign)
library(tibble)




#### 1. Load low throughput (LTP) PPI data from BioGrid and get the unique genes/proteins

## 1.a Download and process BioGRID data
# Get the URL and store it in a temp file
biogrid_url <- "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.244/BIOGRID-ORGANISM-4.4.244.tab3.zip"
temp_zip <- tempfile(fileext = ".zip")

# Download the file using binary mode
options(timeout = 600)  
download.file(biogrid_url, temp_zip, mode = "wb")

# Create a vector with all E. coli strains to be processed
ecoli_strains <- c(
  "Escherichia_coli_K12",
  "Escherichia_coli_K12_MC4100_BW2952",
  "Escherichia_coli_K12_MG1655",
  "Escherichia_coli_K12_W3110"
)

# Function to read and process the data for each E. coli strain 
read_biogrid_data <- function(strain) {
  filename <- paste0("BIOGRID-ORGANISM-", strain, "-4.4.244.tab3.txt")
  unzip(temp_zip, files = filename, exdir = tempdir())
  read.delim(file.path(tempdir(), filename), header = TRUE, stringsAsFactors = FALSE)
}

## 1.b Combine all interactions into one merged dataset
merged_biogrid_data <- map_dfr(ecoli_strains, read_biogrid_data)

## 1.c Filter for LT interactions (positive gold standards)
positive_gold_standards <- merged_biogrid_data %>%
  filter(Throughput == "Low Throughput",
         str_starts(Organism.Name.Interactor.A, "Escherichia coli"),
         str_starts(Organism.Name.Interactor.B, "Escherichia coli")) %>%
  dplyr::select(Official.Symbol.Interactor.A, 
                Official.Symbol.Interactor.B,
                SWISS.PROT.Accessions.Interactor.A,
                SWISS.PROT.Accessions.Interactor.B,
                Entrez.Gene.Interactor.A,
                Entrez.Gene.Interactor.B,
                Organism.ID.Interactor.A,
                Organism.ID.Interactor.B) %>%
  dplyr::rename(Interactor_A = Official.Symbol.Interactor.A, 
                Interactor_B = Official.Symbol.Interactor.B,
                Uniprot_ID_A = SWISS.PROT.Accessions.Interactor.A,
                Uniprot_ID_B = SWISS.PROT.Accessions.Interactor.B,
                Entrez_Gene_ID_A = Entrez.Gene.Interactor.A,
                Entrez_Gene_ID_B = Entrez.Gene.Interactor.B,
                Organism_ID_A = Organism.ID.Interactor.A,
                Organism_ID_B = Organism.ID.Interactor.B)

## 1.d Extract unique proteins/genes
positive_unique_genes <- positive_gold_standards %>%
  dplyr::select(Interactor_A, Uniprot_ID_A, Entrez_Gene_ID_A, Organism_ID_A) %>%
  dplyr::rename(Gene = Interactor_A, 
                Uniprot_ID = Uniprot_ID_A,
                Entrez_Gene_ID = Entrez_Gene_ID_A,
                Organism_ID = Organism_ID_A) %>%
  bind_rows(
    positive_gold_standards %>%
      dplyr::select(Interactor_B, Uniprot_ID_B, Entrez_Gene_ID_B, Organism_ID_B) %>%
      dplyr::rename(Gene = Interactor_B, 
                    Uniprot_ID = Uniprot_ID_B,
                    Entrez_Gene_ID = Entrez_Gene_ID_B,
                    Organism_ID = Organism_ID_B)
  ) %>%
  distinct(Gene, .keep_all = TRUE)

## 1.e Clean temporary files
unlink(temp_zip)
unlink(file.path(tempdir(), "BIOGRID-ORGANISM-*.txt"))





#### 2. Load the negative gold standards and get the unique proteins

## 2.1 Read the csv file (ic_3-5)
ng_file <- read.csv("C:/Users/User/Desktop/Protein interactions/Data collection/Get common proteins - positive, negative, datasets/stella_ecoli_ngs/stella_ecoli_ngs/ngsr_ic_3-5.csv")

# Rename the list and the columns
ng_ic_3_5_data <- ng_file %>%
  dplyr::rename(Interactor_A = Gene.A, Interactor_B = Gene.B)

## 2.2 Get the unique genes from the ng_ic_3_5_data
ng_unique_genes <- ng_ic_3_5_data %>%
  dplyr::select(Interactor_A) %>%
  dplyr::rename(Gene = Interactor_A) %>%
  bind_rows(
    ng_ic_3_5_data %>%
      dplyr::select(Interactor_B) %>%
      dplyr::rename(Gene = Interactor_B)
  ) %>%
  distinct(Gene)

## 2.3 Load and process Ecoli_genes data from json format
file_path <- "C:\\Users\\User\\Desktop\\Protein interactions\\Data collection\\Get common proteins - positive, negative, datasets\\stella_ecoli_ngs\\stella_ecoli_ngs\\ecoli_genes.json"

ng_ecoli_genes <- fromJSON(file_path) %>%
  {data.frame(Gene = names(.), Uniprot_ID = unlist(.), stringsAsFactors = FALSE)}


## 2.4 Find overlaps with E. coli genes
common_between_ng_and_reference <- inner_join(ng_unique_genes, 
                                              ng_ecoli_genes, by = "Gene")

## 2.5 Add the IDs in the  ng interactions table
ng_ic_3_5_with_IDs <- ng_ic_3_5_data %>%
  left_join(ng_ecoli_genes, by = c("Interactor_A" = "Gene")) %>%
  dplyr::rename(Uniprot_ID_A = Uniprot_ID) %>%  
  left_join(ng_ecoli_genes, by = c("Interactor_B" = "Gene")) %>%
  dplyr::rename(Uniprot_ID_B = Uniprot_ID)  





#### 3. Get the list of common proteins between positive and negative gold standards
intersections_ps_ng <- positive_unique_genes %>%
  semi_join(common_between_ng_and_reference, by = "Gene") %>%
  dplyr::select(Gene) 




#### 4. Get the tables of positive and negative interactions that contain the common genes
# Function to filter rows where values of columns A and B are found in the unique genes
filter_interactions <- function(data, col_a, col_b, genes) {
  data %>%
    filter({{col_a}} %in% genes | {{col_b}} %in% genes)
}

## Filter the NG and PS sets for interactions of intersecting genes                                                        
ps_ng_intersect_interactions <- list(
  ps = positive_gold_standards %>% 
    filter_interactions(Interactor_A, Interactor_B, intersections_ps_ng$Gene),
  ng =  ng_ic_3_5_with_IDs%>% 
    filter_interactions(Interactor_A, Interactor_B, intersections_ps_ng$Gene))





#### 5 Load the HTP datasets and find unique proteins of each dataset

## 5.1 Define the path and create the output folder
zip_folder <- "C:/Users/User/Desktop/Protein interactions/Data collection/Get common proteins - positive, negative, datasets"
output_folder <- file.path(zip_folder, "Unzipped_Files")
dir.create(output_folder, showWarnings = FALSE) # prevent warnings of folder already exists

# Define file to unzip
files_to_unzip <- c("BIOGRID-PUBLICATION-182729-4.4.244.DOWNLOADS.zip",
                    "BIOGRID-PUBLICATION-206197-4.4.244.DOWNLOADS.zip",
                    "BIOGRID-PUBLICATION-248142-4.4.244.DOWNLOADS.zip")

# Unzip the files, walk function is applied to each of the files
walk(file.path(zip_folder, files_to_unzip), ~unzip(.x, exdir = output_folder))

# Load all tab files
data_files <- list.files(output_folder, pattern = "\\.(txt|tab)$", full.names = TRUE)
walk(data_files, ~assign(tools::file_path_sans_ext(basename(.x)), fread(.x), envir = .GlobalEnv))


## 5.2 Filter the datasets
# Create a list with filtered datasets
filtered_biog_datasets <- list(
  rajagopala = `BIOGRID-PUBLICATION-182729-4.4.244.tab3`,
  babu = `BIOGRID-PUBLICATION-206197-4.4.244.tab3`,
  arifuzzaman = `BIOGRID-PUBLICATION-248142-4.4.244.tab3`
) %>%
  # Apply filtering to all datasets
  map(~ .x %>%
        filter(Throughput == "High Throughput",
               str_starts(`Organism Name Interactor A`, "Escherichia coli"),
               str_starts(`Organism Name Interactor B`, "Escherichia coli")) %>%
        dplyr::rename(Interactor_A = `Official Symbol Interactor A`, 
                      Interactor_B = `Official Symbol Interactor B`,
                      Uniprot_ID_A = `SWISS-PROT Accessions Interactor A`,
                      Uniprot_ID_B = `SWISS-PROT Accessions Interactor B`,
                      Entrez_Gene_ID_A =`Entrez Gene Interactor A`,
                      Entrez_Gene_ID_B = `Entrez Gene Interactor B`,
                      Organism_ID_A = `Organism ID Interactor A`,
                      Organism_ID_B = `Organism ID Interactor B`)%>%
        dplyr::select(Interactor_A, 
                      Interactor_B, 
                      Uniprot_ID_A, 
                      Uniprot_ID_B,
                      Entrez_Gene_ID_A,
                      Entrez_Gene_ID_B,
                      Organism_ID_A,
                      Organism_ID_B)
  )



##  5.3 Get unique HTP interactions for each dataset
unique_htp_interactions_list <- lapply(filtered_biog_datasets, function(df) {
  # Convert to data.table for speed (if not already)
  dt <- as.data.table(df)
  
  # Create standardised pair with alphabetical order
  dt[, sorted_pair := paste(fifelse(Interactor_A < Interactor_B, 
                                    Interactor_A, Interactor_B),
                            fifelse(Interactor_A < Interactor_B, 
                                    Interactor_B, Interactor_A),
                            sep = "|")]
  
  # Keep first occurrence of each unique pair
  unique_dt <- dt[, .SD[1], by = sorted_pair]
  
  # Remove temporary column and return
  unique_dt[, sorted_pair := NULL]
  return(unique_dt)
})




## 5.4 Get the unique genes from Biogrid data
# Function to extract unique genes from the filtered datasets
extract_unique_genes <- function(df) {
  # Convert to data frame 
  df <- as.data.frame(df)
  
  # Extract genes and IDs for interactor A
  df_a <- data.frame(
    Gene = df[, 1],        # First column (Gene A)
    Uniprot_ID = df[, 3], 
    Entrez_Gene_ID = df[, 5],     
    Organism_ID = df[, 7],
    stringsAsFactors = FALSE
  )
  
  # Extract genes and IDs for interactor B  
  df_b <- data.frame(
    Gene = df[, 2],        # Second column (Gene B)
    Uniprot_ID = df[, 4],  
    Entrez_Gene_ID = df[, 6],     
    Organism_ID = df[, 8],
    stringsAsFactors = FALSE
  )
  
  # Combine and get the unique genes
  unique_genes <- rbind(df_a, df_b) %>% 
    dplyr::distinct(Gene, .keep_all = TRUE) %>%
    dplyr::arrange(Gene)
  
  return(unique_genes)
}



# Apply the function to all filtered biogrid datasets to get the unique genes with IDs
unique_genes_biogrid_datasets <- purrr::map(
  filtered_biog_datasets, 
  extract_unique_genes
) %>%
  set_names(paste0(names(filtered_biog_datasets), "_unique_genes"))






#### 6 Find intersections between positive, negative and each dataset and then intersections across all datasets

# Find intersections with positive-negative sets, for each of the datasets
# imap function loops over in the protein_sets list
intersections_with_ps_ng <- imap(unique_genes_biogrid_datasets, ~ {
  # Keep only rows where Gene matches intersections_ps_ng$Gene
  matched_genes <- .x %>% 
    filter(Gene %in% intersections_ps_ng$Gene) %>% 
    dplyr::select(Gene, Uniprot_ID, Entrez_Gene_ID, Organism_ID)  # Select the 3 columns
  
  return(matched_genes)
})

# Rename the data frames
names(intersections_with_ps_ng) <- paste0(names(unique_genes_biogrid_datasets), "_intersections_with_ps_ng")



# Find intersections across all datasets and ps_ng
# Extract Gene columns from all datasets
gene_vectors <- c(
  list(intersections_ps_ng$Gene),  # Reference genes (already a vector)
  purrr::map(unique_genes_biogrid_datasets, ~ .x$Gene)  # Extract Gene column from each dataset
)

# Find common genes across all datasets
common_genes <- purrr::reduce(gene_vectors, intersect)

# Find common genes across all vectors
intersections_across_datasets <- tibble(Gene = common_genes)





#### 7 Get interactions involving the intersections across all datasets, ps and ng

# Combine all datasets into a list 
all_datasets <- c(
  filtered_biog_datasets,  
  list(ps_data = positive_gold_standards, 
       ng_data = ng_ic_3_5_with_IDs)  
)

# Function to find interactions that involve the intersections across the dataset, ps and ng
find_intersect_interactions <- function(dataset) {
  dataset %>%
    filter(Interactor_A %in% intersections_across_datasets$Gene | 
             Interactor_B %in% intersections_across_datasets$Gene)
}

# Apply the find_intersect_interactions function
int_interactions_across_datasets <- lapply(all_datasets, find_intersect_interactions)

# Add suffix to data frames
names(int_interactions_across_datasets) <- paste0(names(int_interactions_across_datasets), "_int_interactions")

rm(all_datasets)





#### 8 Find the unique pairs of interactors for each data set
# it treats the pairs A->B and B->A as the same pairs
#if one pair has its reverse in the dataset, it will give the first instance of this pair only once
# Function to standardise interaction pairs in alphabetical order
standardise_interaction <- function(df) {
  df %>%
    mutate(
      # Sort interactor names alphabetically and combine them
      sorted_pair = ifelse(
        tolower(Interactor_A) < tolower(Interactor_B),
        paste(Interactor_A, Interactor_B, sep = "_"),
        paste(Interactor_B, Interactor_A, sep = "_")
      )
    ) %>%
    # Keep only the first instance of each unique pair
    distinct(sorted_pair, .keep_all = TRUE) %>%
    # Remove the temporary column
    dplyr::select(-sorted_pair)
}

# Get unique interacting pairs for each dataset
unique_interactions <- int_interactions_across_datasets %>%
  map(standardise_interaction) %>%
  set_names(c("raj", "bab", "arif", "ps_data", "ng_data"))



# Create a copy of the data before any modification
unique_interactions_original_data <- purrr::map(unique_interactions, ~ .x)





# Maybe it wont be used
# #### Find the overlaps between ps and HTP and then the overlaps between ng and HTP
# # Create a key function
# make_interaction_key <- function(df) {
#   df %>%
#     mutate(
#       interaction_key = if_else(
#         Interactor_A < Interactor_B,
#         paste(Interactor_A, Interactor_B, sep = "|"),
#         paste(Interactor_B, Interactor_A, sep = "|")
#       )
#     )
# }
# 
# # List of your three data frames to compare with df4
# df_list <- list(raj = unique_interactions$raj,
#                 bab = unique_interactions$bab,
#                 arif = unique_interactions$arif
#                 )
# 
# # Standardise all data frames (including df4)
# standardised_dfs <- lapply(df_list, make_interaction_key)
# standardised_ps_df <- make_interaction_key(unique_interactions$ps_data)
# standardised_ng_df <- make_interaction_key(unique_interactions$ng_data)
# 
# 
# # Initialise an empty list to store results
# overlapping_interactions_with_ps <- list()
# 
# # Loop through each data frame and compare with df4
# for (i in seq_along(standardised_dfs)) {
#   current_df <- standardised_dfs[[i]]
#   
#   # Find common interactions (keeping all columns from the current df)
#   common_interactions <- current_df %>%
#     filter(interaction_key %in% standardised_ps_df$interaction_key) %>%
#     dplyr::select(-interaction_key)  # Remove the helper key
#   
#   # Store in the list with a name (e.g., "df1_vs_df4")
#   overlapping_interactions_with_ps[[paste0(names(df_list)[i], "_overlap_with_ps")]] <- common_interactions
# }
# 
# 
# # Get the overlapping list of pairs for HTP and NGS
# # Initialise an empty list to store results
# overlapping_interactions_with_ng <- list()
# 
# # Loop through each data frame and compare with df4
# for (i in seq_along(standardised_dfs)) {
#   current_df <- standardised_dfs[[i]]
#   
#   # Find common interactions (keeping all columns from the current df)
#   common_interactions <- current_df %>%
#     filter(interaction_key %in% standardised_ng_df$interaction_key) %>%
#     dplyr::select(-interaction_key)  # Remove the helper key
#   
#   # Store in the list with a name (e.g., "df1_vs_df4")
#   overlapping_interactions_with_ng[[paste0(names(df_list)[i], "_overlap_with_ng")]] <- common_interactions
# }
# 
# rm(common_interactions)
# 



#### 10 Get the table with the numbers of interacting pairs for each of the 84 genes across datasets

# Rename the unique interaction datasets
names(unique_interactions) <- c(
  "raj",  
  "bab",      
  "arif", 
  "ps_data",     
  "ng_data"      
)


# Get interaction data from the nested list
get_interaction_data <- function(dataset_name) {
  
  # remove unwanted columns
  unique_interactions[[dataset_name]] %>% 
    dplyr::select(-matches("Uniprot_ID_A"),
                  -matches("Uniprot_ID_B"),
                  -matches("Entrez_Gene_ID_A"),
                  -matches("Entrez_Gene_ID_B"),
                  -matches("Organism_ID_A"), 
                  -matches("Organism_Name_Interactor_A"),
                  -matches("Organism_Name_Interactor_B"))
}

# List all dataset names to be analysed
dataset_names <- c("raj",  
                   "bab",      
                   "arif", 
                   "ps_data",     
                   "ng_data") 


# Create the function that creates bidirectional pairs and counts unique partners per gene
process_interactions <- function(interaction_data, dataset_name) {
  intersections_across_datasets %>%
    left_join(
      bind_rows(
        interaction_data %>% transmute(gene = Interactor_A, partner = Interactor_B),
        interaction_data %>% transmute(gene = Interactor_B, partner = Interactor_A)
      ) %>%
        distinct() %>%
        group_by(gene) %>%
        summarize(
          !!paste0("No._partners_", dataset_name) := n(),
          !!paste0("Partners_", dataset_name) := toString(sort(unique(partner))),
          .groups = "drop"),
      by = c("Gene" = "gene")
    ) 
}



# Process all datasets and merge results
combined_results <- purrr::reduce(
  map(dataset_names, ~ {
    interaction_data <- get_interaction_data(.x)
    process_interactions(interaction_data, .x)
  }),
  full_join,
  by = "Gene"
)


# Get total number of partners across all datasets
total_interacting_partners <- combined_results %>%
  mutate(
    Total_partners = rowSums(dplyr::select(., starts_with("No._partners_")))
  ) %>%
  # Select the columns for the table
  dplyr::select(Gene, starts_with("No._partners_"), starts_with("Partners_"), 
                Total_partners)


rm(combined_results)







############### Selection of Interactions to be analysed
#### 12 AlphaFold data selection 
#### Random data selection to be run on AlphaFold

# List all datasets to bi filtered
data_to_be_randomly_selected <- list(unique_interactions$raj, 
                                     unique_interactions$bab, 
                                     unique_interactions$arif,
                                     unique_interactions$ng_data)

# Name the datasets (optional but helpful for tracking)
names(data_to_be_randomly_selected) <- c("raj", "bab", "arif", "ng_data")

# Set sample size of 130 interactions as per dataset
sample_size <- 130  

# Randomly select 130 interactions from each dataset (
set.seed(123) 
random_interactions <- lapply(data_to_be_randomly_selected, function(df) {
  df %>% slice_sample(n = sample_size)  # Randomly selects 130 rows
})

rm(data_to_be_randomly_selected)

# # # write the data in an excel sheet
# library(writexl)
# write_xlsx(positive_gold_standards, "positive_gold_standards.xlsx")
# # 
# 





##### 13 Get the protein sequence from using the API Uniprot address for the randomly selected data and the ps data
# this API searched by gene names (GN) and the taxonomy ID of the organism
## there were issues parsing the sequences for substrains so the sequences will be obtained using
# the "parent" organism ID (83333)

# Get all the data to be run on AlphaFold together
alphafold_list <- list(
  raj = random_interactions$raj,
  bab = random_interactions$bab,
  arif = random_interactions$arif,
  ng_data = random_interactions$ng_data,
  ps_data = unique_interactions$ps_data
)


# Get the unique genes from all sets into a vector
all_genes <- unique(unlist(lapply(alphafold_list, function(df) {
  c(df$Interactor_A, df$Interactor_B)
})))



# Function to get sequence from uniprot API
get_protein_sequence <- function(gene, organism_id = "83333", proteome_id = "UP000000625") {
  library(httr)
  library(Biostrings)
  library(R.utils)
  
  # UniProt API URL with the following parameters
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=(",
    "gene:", gene,
    "+AND+reviewed:true",
    "+AND+organism_id:", organism_id,
    "+AND+proteome:", proteome_id,
    ")"
  )
  
  # Create temp files to store compressed and decompressed FASTA
  gz_file <- tempfile(fileext = ".fasta.gz")
  fasta_file <- tempfile(fileext = ".fasta")
  
  # Download and decompress
  response <- GET(url, write_disk(gz_file, overwrite = TRUE))
  if (response$status_code != 200) {
    return(list(sequence = NA, 
                status = "api_error", 
                length = NA,
                uniprot_id = NA_character_))
  }
  
  gunzip(gz_file, fasta_file, overwrite = TRUE)
  
  # Parse FASTA
  fasta <- readAAStringSet(fasta_file)
  
  # Handle different cases
  if (length(fasta) == 0) {
    return(list(sequence = NA, status = "no_sequence", length = NA))
  } else if (length(fasta) > 1) {
    # Extract all UniProt IDs for multiple sequences
    uniprot_ids <- sapply(names(fasta), function(x) {
      gsub("^.*\\|([A-Z0-9]+)\\|.*$", "\\1", x)
    })
    return(list(sequence = NA, status = "multiple_sequences", length = NA,
                uniprot_id = paste(uniprot_ids, collapse = ",")))
  } else {
    seq_string <- as.character(fasta[[1]])
    # Extract UniProt ID from FASTA header (e.g., ">sp|P0A6F5|DNAA_ECOLI" → "P0A6F5")
    uniprot_id <- gsub("^.*\\|([A-Z0-9]+)\\|.*$", "\\1", names(fasta)[1])
    return(list(sequence = seq_string, status = "success", length = nchar(seq_string), uniprot_id = uniprot_id))
  }
}




# Function to process protein interaction dataframe
process_genes <- function(genes, organism_id = "83333", proteome_id = "UP000000625") {
  library(dplyr)
  
  # Initialize results dataframe with additional status columns
  results <- data.frame(
    gene = genes,
    sequence = character(length(genes)),
    sequence_length = integer(length(genes)),
    uniprot_id = character(length(genes)),
    status = character(length(genes)),
    stringsAsFactors = FALSE
  )
  
  # Process each gene
  for (i in seq_along(genes)) {
    cat(sprintf("Processing gene %d/%d: %s\n", i, length(genes), genes[i]))
    
    result <- tryCatch({
      get_protein_sequence(genes[i], organism_id, proteome_id)
    }, error = function(e) {
      list(sequence = NA, status = "error", length = NA, uniprot_id = NA_character_)
    })
    
    results$sequence[i] <- result$sequence %||% NA_character_
    results$sequence_length[i] <- result$length %||% NA_integer_
    results$uniprot_id[i] <- result$uniprot_id %||% NA_character_
    results$status[i] <- result$status
    Sys.sleep(0.5)
  }
  
  return(results)
}
  




# Function to filter interactions and merge sequences
filter_and_merge_interactions <- function(interaction_df, sequence_df) {
  library(dplyr)
  
  # Get successful genes only
  successful_genes <- sequence_df %>% 
    filter(status == "success") %>% 
    pull(gene)
  
  # Filter interactions to only include successful genes
  filtered_interactions <- interaction_df %>%
    filter(
      Interactor_A %in% successful_genes & 
        Interactor_B %in% successful_genes
    )
  
  # Merge sequences
  final_df <- filtered_interactions %>%
    left_join(
      sequence_df %>% 
        dplyr::select(gene, uniprot_id, sequence, sequence_length) %>%
        dplyr::rename(Sequence_A = sequence, 
                      Uniprot_ID_seq_A = uniprot_id,
                      Length_A = sequence_length),
      by = c("Interactor_A" = "gene")
    ) %>%
    left_join(
      sequence_df %>% 
        dplyr::select(gene, uniprot_id, sequence, sequence_length) %>%
        dplyr::rename(Sequence_B = sequence, 
                      Uniprot_ID_seq_B = uniprot_id,
                      Length_B = sequence_length),
      by = c("Interactor_B" = "gene")
    )
  
  # Summary statistics to report the number of interactions removed due to multiple or missing sequences
  original_count <- nrow(interaction_df)
  final_count <- nrow(final_df)
  removed_count <- original_count - final_count
  
  cat(sprintf("Filtered %d interactions (%d interactions removed)\n", 
              final_count, removed_count, (removed_count/original_count)*100))
  
  return(list(
    interactions = final_df,
    stats = list(
      original = original_count,
      final = final_count,
      removed = removed_count,
      removal_rate = (removed_count/original_count)*100
    )
  ))
}



process_interaction_dataframes <- function(df_list, df_names = NULL, organism_id = "83333", proteome_id = "UP000000625") {
  
  # Get all unique genes
  all_genes <- unique(unlist(lapply(df_list, function(df) {
    c(df$Interactor_A, df$Interactor_B)
  })))
  
  cat("=== PROCESSING SUMMARY ===\n")
  cat("Dataframes:", length(df_list), "\n")
  cat("Total interactions:", sum(sapply(df_list, nrow)), "\n")
  cat("Unique genes:", length(all_genes), "\n\n")
  
  # Get sequences for all genes
  sequences <- process_genes(all_genes, organism_id, proteome_id)
  
  # Print sequence summary
  status_summary <- table(sequences$status)
  cat("=== SEQUENCE RETRIEVAL SUMMARY ===\n")
  print(status_summary)
  cat("Success rate:", round(status_summary["success"]/length(all_genes)*100, 1), "%\n\n")
  
  # Process each dataframe
  results <- list()
  for (i in seq_along(df_list)) {
    cat("Processing", df_names[i], "...\n")
    result <- filter_and_merge_interactions(df_list[[i]], sequences)
    results[[df_names[i]]] <- result$interactions
    cat("\n")
  }
  
  return(list(
    dataframes = results,
    sequences = sequences,
    genes_excluded = sequences %>% filter(status != "success") %>% pull(gene)
  ))
}



df_names <- c("raj", "bab", "arif", "ng_data", "ps_data")


alphafold_results <- process_interaction_dataframes(alphafold_list, df_names)
alphafold_list_with_seq_final <- alphafold_results$dataframes


# write the file
library(writexl)
write_xlsx(alphafold_list_with_seq_final, path = "alphafold_list_with_seq_final.xlsx")



## Check the removed pairs from each data set

alphafold_removed_pairs_list <- lapply(names(alphafold_list), function(name) {
  original_df <- alphafold_list[[name]]
  filtered_df <- alphafold_list_with_seq_final[[name]]
  
  # Identify removed pairs
  removed_pairs <- anti_join(
    original_df,
    filtered_df,
    by = c("Interactor_A", "Interactor_B")
  )
  
  # Add dataframe name for reference
  removed_pairs$source_df <- name
  return(removed_pairs)
})

# Combine all removed pairs into one dataframe
all_removed_pairs <- bind_rows(alphafold_removed_pairs_list)

# Print summary
cat("=== Total Removed Interactions Across All Dataframes ===\n")
print(nrow(all_removed_pairs))

cat("\n=== Removed Pairs by Dataframe ===\n")
print(all_removed_pairs %>% count(source_df))

# View detailed removed pairs
cat("\n=== Full List of Removed Pairs ===\n")
print(all_removed_pairs)








####### 14 Get the FASTA files to run Hmmer
#### 14 Prepare the data to run on Hmmer using FASTA sequence

## Get sequences for yajL 57 interacting partners found in PS data
# 1 Add source ID for each of the data frames htp, ng and ps
unique_interactions$raj$Source <- "raj"
unique_interactions$bab$Source <- "bab"
unique_interactions$arif$Source <- "arif"
unique_interactions$ps_data$Source <- "ps_data"
unique_interactions$ng_data$Source <- "ng_data"

# 2 Put all data and sources together
all_sources_combined <- bind_rows(unique_interactions$raj,
                                  unique_interactions$bab,
                                  unique_interactions$arif,
                                  unique_interactions$ng_data,
                                  unique_interactions$ps_data)



## 3 Get the interactors of yajL from PS sets
yajL_ps_parnters <- unique_interactions$ps_data %>% 
  filter(Interactor_A == "yajL" | Interactor_B == "yajL") %>%
  mutate(Interactor = if_else(Interactor_A == "yajL", Interactor_B, Interactor_A),
         Uniprot_ID = if_else(Interactor_A == "yajL", Uniprot_ID_B, Uniprot_ID_A),
         Entrez_Gene_ID = if_else(Interactor_A == "yajL", Entrez_Gene_ID_B, Entrez_Gene_ID_A)) %>%
  dplyr::select(Interactor, 
                Uniprot_ID,
                Entrez_Gene_ID,
                Source)

# Identify the unique ps interacting partners
yajL_unique_ps_parnters <- unique(yajL_ps_parnters$Interactor) 



## 4 Find if any of the 57 yajL interacting partners  are found to interact with yajL in any other data sets
# Get all YajL interactors across ALL datasets (with sources)
yajL_all_interactors <- all_sources_combined %>%
  filter(Interactor_A == "yajL" | Interactor_B == "yajL") %>%
  mutate(
    Interactor = if_else(Interactor_A == "yajL", Interactor_B, Interactor_A),
    Uniprot_ID = if_else(Interactor_A == "yajL", Uniprot_ID_B, Uniprot_ID_A),
    Entrez_Gene_ID = if_else(Interactor_A == "yajL", Entrez_Gene_ID_B, Entrez_Gene_ID_A)
  ) %>%
  dplyr::select(Interactor, Uniprot_ID, Entrez_Gene_ID, Source)

# Filter for the 57 PS partners and find their presence in other data sets
yajL_ps_partners_all_sources <- yajL_all_interactors %>%
  # Keep only the 57 PS partners (assuming yajL_unique_ps_parnters contains their names)
  filter(Interactor %in% yajL_unique_ps_parnters) %>%
  # Group by interactor to collapse sources
  group_by(Interactor, Uniprot_ID, Entrez_Gene_ID) %>%
  summarise(
    Source = toString(sort(unique(Source))),  # Combine all sources into one string
    .groups = "drop"
  )


# Get unique proteins but multiple sources if exist
yajL_57_partners_all_sources_combined <- yajL_ps_partners_all_sources %>%
  group_by(Interactor) %>%
  # combine data together separated by "," for each column
  dplyr::summarize(
    Source = paste(unique(Source), collapse = ", "),
    .groups = "drop"
  ) %>%
  ungroup() 

yajL_57_partners_all_sources_combined <- yajL_57_partners_all_sources_combined %>%
  dplyr::select(-Uniprot_ID, -Entrez_Gene_ID)




#### get the Fasta sequences for yajL 57 ps interactors
# Function to process interactors and save sequences with custom headers
process_57_ps_interactors <- function(interactor_df, 
                                organism_id = "83333",
                                proteome_id = "UP000000625",
                                output_filename = "yajL_57_ps_interactor_sequences.fasta") {
  
  # Load required packages
  library(dplyr)
  library(purrr)
  
  # Initialize results list
  results <- list(
    fasta_entries = character(),
    failed_entries = tibble(
      Interactor = character(),
      Reason = character()
    )
  )
  
  # Process each interactor
  walk(1:nrow(interactor_df), function(i) {
    gene <- interactor_df$Interactor[i]
    source <- interactor_df$Source[i]
    
    # Get sequence using existing function
    seq_result <- get_protein_sequence(gene, organism_id, proteome_id)
    
    if (seq_result$status == "success") {
      # Create custom FASTA header
      header <- sprintf(">%s | %s | %s",
                        gene,
                        seq_result$uniprot_id,  # Use the UniProt ID from the fetched sequence
                        source)
      
      # Add to results
      results$fasta_entries <<- c(results$fasta_entries,
                                  paste0(header, "\n", seq_result$sequence))
    } else {
      # Record failure
      results$failed_entries <<- add_row(results$failed_entries,
                                         Interactor = gene,
                                         Reason = seq_result$status)
    }
  })
  
  # Save successful sequences to file
  if (length(results$fasta_entries) > 0) {
    writeLines(results$fasta_entries, output_filename)
    message(sprintf("Saved %d sequences to %s", 
                    length(results$fasta_entries), 
                    output_filename))
  } else {
    warning("No sequences were successfully retrieved")
  }
  
  # Print summary
  if (nrow(results$failed_entries) > 0) {
    message("\nFailed retrievals:")
    print(count(results$failed_entries, Reason, name = "Count"))
  }
  
  # Return invisible list with all data
  invisible(list(
    fasta_file = if (length(results$fasta_entries) > 0) output_filename else NULL,
    successful_sequences = results$fasta_entries,
    failed_retrievals = results$failed_entries
  ))
}

# Save the results for 57 yajL interacting partners
yajL_57_ps_interactors_fasta_seq_resutls <- process_57_ps_interactors(yajL_57_partners_all_sources_combined)






### Get the fasta files of sequences for the interactors that are included in HTP and the PS data sets

get_gene_interactors_and_fasta <- function(gene_name, 
                                           unique_interactions,
                                           organism_id = "83333",
                                           proteome_id = "UP000000625",
                                           output_filename = NULL) {
  
  # Load required packages
  library(dplyr)
  library(purrr)
  
  # Set default output filename if not provided
  if (is.null(output_filename)) {
    output_filename <- paste0(gene_name, "_htp_ps_interactors.fasta")
  }
  
  cat(sprintf("Analysing interactors for gene: %s\n", gene_name))
  
  ## Helper function to create FASTA header
  create_fasta_header <- function(interactor, uniprot_id, gene_name, source) {
    sprintf(">%s | %s | %s_interactor | %s", 
             interactor, uniprot_id, gene_name, source)
  }
  
  ## Helper function to extract interactors from a dataset
  extract_interactors <- function(df, gene_name, source_name) {
    if (inherits(df, "data.table")) data.table::setDF(df)
    
    df %>%
      filter(Interactor_A == gene_name | Interactor_B == gene_name) %>%
      mutate(
        Interactor = if_else(Interactor_A == gene_name, Interactor_B, Interactor_A),
        Source = source_name
      ) %>%
      dplyr::select(Interactor, Source)
  }
  
  ## 1. Extract interactors from all datasets
  # Process HT and PS datasets
  ht_ps_interactors <- list(
    raj = extract_interactors(unique_interactions$raj, gene_name, "raj"),
    bab = extract_interactors(unique_interactions$bab, gene_name, "bab"),
    arif = extract_interactors(unique_interactions$arif, gene_name, "arif"),
    ps_data = extract_interactors(unique_interactions$ps_data, gene_name, "ps_data")
  ) %>%
    bind_rows()
  
  # Process NG dataset to get interactor list
  ng_interactors <- unique_interactions$ng_data %>%
    filter(Interactor_A == gene_name | Interactor_B == gene_name) %>%
    mutate(Interactor = if_else(Interactor_A == gene_name, Interactor_B, Interactor_A)) %>%
    distinct(Interactor) %>%
    pull(Interactor)
  
  # Check if any interactors were found
  if (nrow(ht_ps_interactors) == 0) {
    cat(sprintf("No interactors found for gene %s in PS/HTP datasets\n", gene_name))
    return(NULL)
  }
  
  ## 2. Combine and consolidate all interactor data
  all_interactors <- ht_ps_interactors %>%
    group_by(Interactor) %>%
    summarize(
      Source = paste(unique(Source), collapse = ", "),
      .groups = "drop"
    ) %>%
    # Add NG data source if interactor is found in NG dataset
    mutate(
      Source = if_else(
        Interactor %in% ng_interactors,
        paste(Source, "ng_data", sep = ", "),
        Source
      )
    )
  
  cat(sprintf("Found %d unique interactors for %s\n", nrow(all_interactors), gene_name))
  
  ## 3. Fetch sequences and process results
  cat("Fetching FASTA sequences using the Uniprot API address...\n")
  
  # Process all interactors and collect results
  sequence_results <- all_interactors %>%
    rowwise() %>%
    mutate(
      seq_result = list(get_protein_sequence(Interactor, organism_id, proteome_id))
    ) %>%
    ungroup() %>%
    mutate(
      status = map_chr(seq_result, ~.x$status),
      sequence = map_chr(seq_result, ~if(.x$status == "success") .x$sequence else NA_character_),
      fetched_uniprot_id = map_chr(seq_result, ~if(.x$status == "success") .x$uniprot_id else NA_character_)) %>%
    dplyr::select(-seq_result)
  
  # Separate valid and excluded interactors
  valid_interactors <- sequence_results %>%
    filter(status == "success") %>%
    mutate(
      fasta_header = create_fasta_header(fetched_uniprot_id, Interactor, gene_name, Source),
      fasta_entry = paste0(fasta_header, "\n", sequence)
    )
  
  excluded_interactors <- sequence_results %>%
    filter(status != "success")
  
  ## 4. Generate output
  if (nrow(valid_interactors) == 0) {
    cat("No valid interactors with unambiguous sequences found.\n")
    return(NULL)
  }
  
  # Write FASTA file
  fasta_content <- paste(valid_interactors$fasta_entry, collapse = "\n")
  writeLines(fasta_content, output_filename)
  
  # Prepare return object
  result <- list(
    fasta_file = output_filename,
    valid_interactors = valid_interactors %>% dplyr::select(-fasta_header, -fasta_entry, -sequence, -fetched_uniprot_id),
    excluded_interactors = excluded_interactors %>% dplyr::select(-sequence, -fetched_uniprot_id),
    exclusion_stats = table(excluded_interactors$status)
  )
  
  # Print summary
  cat(sprintf("\n=== Results Summary ===\n"))
  cat(sprintf("Valid interactors with sequences: %d\n", nrow(valid_interactors)))
  cat(sprintf("Excluded interactors: %d\n", nrow(excluded_interactors)))
  if (nrow(excluded_interactors) > 0) {
    cat("Exclusion reasons:\n")
    print(result$exclusion_stats)
  }
  cat(sprintf("Results saved to: %s\n", output_filename))
  
  return(result)
}


# Greate a list of genes I want the fasta files of their HTP and PS interactors
gene_list <- c("yajL", 
               "clpB", 
               "crp", 
               "dnaX",
               "minC",
               "clpA",
               "dnaJ",
               "dnaN",
               "fldA",
               "frr")


# Create a function to loop over the genes and get the fasta files of their interactors
get_multiple_interactor_fasta_files <- function(gene_list, unique_interactions) {
  results <- list()
  
  for (gene in gene_list) {
    cat(sprintf("\n=== Processing gene: %s ===\n", gene))
    
    # Store the result for each gene
    result <- get_gene_interactors_and_fasta(gene, unique_interactions)
    
    # Add to results list (even if NULL)
    results[[gene]] <- result
  }
  
  return(results)
}

htp_ps_partners_fasta_results <- get_multiple_interactor_fasta_files(gene_list, unique_interactions)




#### Fetch the sequences of these 10 genes and save it in Fasta files
yajL_fasta_file <- get_protein_sequence("yajL")


# Function to fetch and save individual FASTA files for a list of genes
fetch_individual_fasta_files <- function(gene_list, 
                                         organism_id = "83333", 
                                         proteome_id = "UP000000625",
                                         output_dir = "intersecting_genes_fasta_files") {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(sprintf("Created directory: %s\n", output_dir))
  }
  
  # Initialize results tracking
  results <- list(
    successful = character(),
    failed = data.frame(
      gene = character(),
      reason = character(),
      uniprot_id = character(),
      stringsAsFactors = FALSE
    )
  )
  
  cat(sprintf("Processing %d genes...\n", length(gene_list)))
  
  # Process each gene
  for (gene in gene_list) {
    cat(sprintf("Fetching sequence for gene: %s\n", gene))
    
    # UniProt API query URL (same as in get_protein_sequence)
    url <- paste0(
      "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=(",
      "gene:", gene,
      "+AND+reviewed:true",
      "+AND+organism_id:", organism_id,
      "+AND+proteome:", proteome_id,
      ")"
    )
    
    # Temp files for download
    gz_file <- tempfile(fileext = ".fasta.gz")
    fasta_file <- tempfile(fileext = ".fasta")
    
    # Download FASTA
    response <- GET(url, write_disk(gz_file, overwrite = TRUE))
    
    if (response$status_code != 200) {
      results$failed <- rbind(results$failed, 
                              data.frame(gene = gene, reason = "api_error"))
      cat(sprintf("  ✗ Failed: API request failed (HTTP %d)\n", response$status_code))
      next
    }
    
    # Decompress and read raw FASTA (before parsing)
    gunzip(gz_file, fasta_file, overwrite = TRUE)
    fasta_lines <- readLines(fasta_file)
    
    # Check if empty
    if (length(fasta_lines) == 0) {
      results$failed <- rbind(results$failed, 
                              data.frame(gene = gene, reason = "no_sequence"))
      cat("  Failed: No sequence found\n")
      next
    }
    
    # Extract header (first line) and sequence (remaining lines)
    fasta_header <- fasta_lines[1]  # Original UniProt header
    fasta_seq <- paste(fasta_lines[-1], collapse = "\n")
    
    # Save to file
    output_file <- file.path(output_dir, paste0(gene, ".fasta"))
    writeLines(c(fasta_header, fasta_seq), output_file)
    
    # Track success
    results$successful <- c(results$successful, output_file)
    cat(sprintf("  ✓ Saved: %s\n", output_file))
  }
  
  # Summary report
  cat("\n=== Summary ===\n")
  cat(sprintf("Successfully processed: %d genes\n", length(results$successful)))
  cat(sprintf("Failed: %d genes\n", nrow(results$failed)))
  
  if (nrow(results$failed) > 0) {
    cat("\nFailure reasons:\n")
    print(table(results$failed$reason))
  }
  
  return(results)
}



test_fasta_files <- fetch_individual_fasta_files(gene_list = gene_list,
                                                 organism_id = "83333",       
                                                 proteome_id = "UP000000625", 
                                                 output_dir = "intersecting_genes_fasta_files" 
)







####### 16 Hmmer output file analysis Domain similarity
##### Hmmer output data analysis for those 57 yajL interacting partners

#### 1 Function to read and parse domtblout files in R
read_domtblout <- function(file) {
  # Read all lines
  lines <- readLines(file)
  # Remove comment lines (starting with #) or empty lines
  data_lines <- lines[!grepl("^#", lines) & lines != ""]
  
  if (length(data_lines) == 0) {
    warning("No data lines found in the file")
    return(data.frame())
  }
  
  # Define column names for domtblout format
  col_names <- c(
    "target_name", "target_accession", "tlen", "query_name", 
    "query_accession", "qlen", "E_value", "score", "bias",
    "domain_num", "domain_of", "c_Evalue", "i_Evalue", 
    "domain_score", "domain_bias", "hmm_from", "hmm_to",
    "ali_from", "ali_to", "env_from", "env_to", "acc", "description"
  )
  
  # Parse each line
  parsed_data <- lapply(data_lines, function(line) {
    # Split by whitespace, but be careful with the description field
    fields <- unlist(strsplit(line, "\\s+"))
    
    # The description field (last field) may contain spaces
    # So we need to handle cases where there are more than 23 fields
    if (length(fields) >= 23) {
      # Combine all fields from position 23 onwards into description
      if (length(fields) > 23) {
        fields[23] <- paste(fields[23:length(fields)], collapse = " ")
      }
      fields <- fields[1:23]
    } else if (length(fields) < 23) {
      # If there are fewer fields, pad with NA
      fields <- c(fields, rep(NA, 23 - length(fields)))
    }
    
    return(fields)
  })
  
  
  # Convert to data frame
  domtbl_df <- do.call(rbind, parsed_data)
  domtbl_df <- as.data.frame(domtbl_df, stringsAsFactors = FALSE)
  colnames(domtbl_df) <- col_names
  
  # Convert numeric columns to appropriate types
  numeric_cols <- c("tlen", "qlen", "E_value", "score", "bias", 
                    "domain_num", "domain_of", "c_Evalue", "i_Evalue", 
                    "domain_score", "domain_bias", "hmm_from", "hmm_to", 
                    "ali_from", "ali_to", "env_from", "env_to", "acc")
  
  for (col in numeric_cols) {
    domtbl_df[[col]] <- as.numeric(domtbl_df[[col]])
  }
  
  return(domtbl_df)
}



# 2 Load HMMER domtblout files
hmmer_57_yajl_ps_prntrs_domtblout <- read_domtblout("C:/Users/User/Desktop/Protein interactions/Data collection/Get common proteins - positive, negative, datasets/hmmscan_results/curated/yajL_57_ps_partners/YAJL_57_ps_partners_curated_pfam_results.domtblout")                 


# Function to read all domtblout files from a directory
read_domtblout_directory <- function(directory_path, pattern = "\\.domtblout$") {
  # Get all domtblout files in the directory
  file_paths <- list.files(directory_path, pattern = pattern, full.names = TRUE)
  
  if (length(file_paths) == 0) {
    warning(paste("No files matching pattern", pattern, "found in", directory_path))
    return(list())
  }
  
  # Read all files
  domtbl_list <- lapply(file_paths, read_domtblout)
  
  # Use basename as names (remove .domtblout extension for cleaner names)
  names(domtbl_list) <- tools::file_path_sans_ext(basename(file_paths))
  
  return(domtbl_list)
}


# Read all domtblout files from intersections directory
domtblout_list_htp_ps_intersections <- read_domtblout_directory(
  "C:/Users/User/Desktop/Protein interactions/Data collection/Get common proteins - positive, negative, datasets/hmmscan_results/curated/intersections/"
)

# Read all domtblout files from htp_ps_partners directory
domtblout_list_htp_ps_prtnts <- read_domtblout_directory(
  "C:/Users/User/Desktop/Protein interactions/Data collection/Get common proteins - positive, negative, datasets/hmmscan_results/curated/htp_ps_partners/"
)





############################################
#### DO the downstream analysis
######### Do domain analysis of interacting partners
## Start with yajL interacting proteins
# 2 Function to get a summary of domain hits for each protein
get_domain_summary <- function(domains_df) {
  domain_summary <- domains_df %>%
    dplyr::select(query_name, target_accession) %>%
    distinct() %>%
    group_by(query_name) %>%
    summarise(domains = paste(unique(target_accession), collapse = ", "),
              num_domains = n_distinct(target_accession),
              .groups = "drop")
  
  return(domain_summary)
}

domain_summary_yajL_57_ps_prtnrs <- get_domain_summary(hmmer_57_yajl_ps_prntrs_domtblout)


# 3 Function that counts domain occurrences across proteins
get_domain_frequency <- function(domains_df) {
  domain_freq <- domains_df %>%
    group_by(target_accession) %>%
    summarise(
      protein_count = n_distinct(query_name),  # Unique proteins with this domain
      total_occurrences = n()                # Total times domain appears (including duplicates)
    ) %>%
    arrange(desc(protein_count))
  
  return(domain_freq)
}

domain_freq_yajL_57_ps_prtnrs <- get_domain_frequency(hmmer_57_yajl_ps_prntrs_domtblout)


# 4 Function to get domains the appear in 2 or more different proteins
filter_shared_domains <- function(domains_df, min_proteins = 2) {
  # Calculate domain frequency
  domain_freq <- get_domain_frequency(domains_df)
  
  # Get domains found in 2 or more proteins
  common_domains <- domain_freq %>%
    filter(protein_count >= min_proteins) %>%
    pull(target_accession)
  
  # Filter data to only include the common domains
  filtered_domains <- domains_df %>%
    filter(target_accession %in% common_domains)
  
  return(list(
    filtered_data = filtered_domains,
    common_domains = common_domains,
    domain_frequency = domain_freq
  ))
}

shared_domains_yajL_57_ps_prtnrs <- filter_shared_domains(hmmer_57_yajl_ps_prntrs_domtblout)


# 5 Function to generate a domain-protein matrix
get_domain_matrix <- function(domains_df, min_proteins = 2) { # the min_proteins = 2 is not used)
  
  # Create matrix
  domain_matrix <- domains_df %>%  
    dplyr::select(query_name, target_accession) %>%
    distinct() %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = target_accession, values_from = value, values_fill = 0) %>%
    as.data.frame()
  
  # Set row names
  rownames(domain_matrix) <- domain_matrix$query_name
  domain_matrix$query_name <- NULL
  
  return(domain_matrix)
}

domain_matrix_shared_proteins_yajL_57_ps_prtnrs <- get_domain_matrix(hmmer_57_yajl_ps_prntrs_domtblout)


# 6 Function to calculate the domain copy number within the same protein
get_domain_copies <- function(domains_df) {
  domain_copies <- domains_df %>%
    group_by(domain = target_accession, protein = query_name) %>%
    summarise(copy_count = n(), .groups = "drop")
  
  # Show proteins with multiple copies of same domain
  multi_copy_domains <- domain_copies %>%
    filter(copy_count > 1) %>%
    arrange(desc(copy_count))
  
  return(list(
    all_copies = domain_copies,
    multi_copy = multi_copy_domains
  ))
}

domain_copies_yajL_57_ps_prtnrs <- get_domain_copies(hmmer_57_yajl_ps_prntrs_domtblout)





### get the heatmaps for yaiL_ps_partners
# heatmap for all domains

# get the matrix
domain_matrix_shared_proteins_yajL_57_ps_prtnrs <- get_domain_matrix(hmmer_57_yajl_ps_prntrs_domtblout)

# Extract headers from the FASTA sequences
yajL_57_ps_fasta_headers <- yajL_57_ps_interactors_fasta_seq_resutls$successful_sequences
# Extract just the header line (everything before the first \n)
header_lines <- str_extract(yajL_57_ps_fasta_headers, "^>[^\\n]+")

# Parse the headers
yajL_57_ps_prtnrs_headers <- str_split_fixed(header_lines, " \\| ", n = 3) %>%
  as.data.frame() %>%
  setNames(c("Protein", "UniProt_ID", "Source")) %>%
  mutate(
    Protein = str_remove(Protein, "^>"),  # Remove the ">" from protein names
    across(everything(), str_trim)
  )

# Create full header string and set as row names
full_headers <- paste(yajL_57_ps_prtnrs_headers$Protein, 
                      yajL_57_ps_prtnrs_headers$UniProt_ID, 
                      yajL_57_ps_prtnrs_headers$Source, 
                      sep = " | ")

# Create a mapping from protein names to full headers
header_mapping <- setNames(full_headers, yajL_57_ps_prtnrs_headers$Protein)

# Map the current row names to full headers
current_rownames <- rownames(domain_matrix_shared_proteins_yajL_57_ps_prtnrs)
rownames(domain_matrix_shared_proteins_yajL_57_ps_prtnrs) <- header_mapping[current_rownames]

# TRANSPOSE the matrix to flip axes (domains on y-axis, proteins on x-axis)
domain_matrix_transposed <- t(domain_matrix_shared_proteins_yajL_57_ps_prtnrs)

# Plot all domains
library(ComplexHeatmap)
ht_yajL_57_ps <- Heatmap(
  domain_matrix_transposed,
  col = c("0" = "grey", "1" = "red"),
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",        
  column_names_side = "bottom",
  column_title = "Domains found in 57 PS interacting partners of yajL",
  row_names_gp = gpar(fontsize = 6),      # Domain label size
  column_names_gp = gpar(fontsize = 8),   # Protein header size
  column_names_rot = 90,                  # Rotated protein names
  heatmap_legend_param = list(
    at = c(0, 1),
    labels = c("Absent", "Present"),
    title = "Domain"
  ),
  cluster_rows = TRUE,           # No row clustering
  cluster_columns = TRUE,        # No column clustering
  show_row_dend = FALSE,          # Explicitly hide row dendrogram
  show_column_dend = FALSE,       # Explicitly hide column dendrogram
  row_order = 1:nrow(domain_matrix_transposed),    # Maintain original row order
  column_order = 1:ncol(domain_matrix_transposed),
  rect_gp = gpar(col = "dimgrey", lwd = 0.5)
)

print(ht_yajL_57_ps)

png("yajL_57_ps_prtnrs_domain_heatmap6.png",
    width = 1600,         # 300 dpi
    height = 2300,
    res = 200)
draw(ht_yajL_57_ps, padding = unit(c(2, 2, 10, 2), "mm"))
dev.off()


#Plot the shared domains only 
# Get shared domains only
shared_domains_yajL_57_ps_prtnrs <- filter_shared_domains(hmmer_57_yajl_ps_prntrs_domtblout)

# Create matrix for only shared domains
domain_matrix_shared_dom_57_ps_yajL_prtnrs <- get_domain_matrix(shared_domains_yajL_57_ps_prtnrs$filtered_data)

# Apply header processing
# keep headers for proteins that actually have shared domains
existing_57_ps_prtnrs <- rownames(domain_matrix_shared_dom_57_ps_yajL_prtnrs) 
header_mapping_shared_dom <- header_mapping[names(header_mapping) %in% existing_57_ps_prtnrs]

# apply the mapping
rownames(domain_matrix_shared_dom_57_ps_yajL_prtnrs) <- header_mapping_shared_dom[rownames(domain_matrix_shared_dom_57_ps_yajL_prtnrs)]

# Transpose the matrix
domain_matrix_shared_dom_57_ps_yajL_prtnrs_transposed <- 
  t(domain_matrix_shared_dom_57_ps_yajL_prtnrs)

# Create heatmap for shared domains
ht_57_ps_shared_domains <- Heatmap(
  domain_matrix_shared_dom_57_ps_yajL_prtnrs_transposed,
  col = c("0" = "grey", "1" = "red"),
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_side = "bottom",
  column_title = "Shared domains in PS interacting partners of yajL ",
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 90,
  heatmap_legend_param = list(
    at = c(0, 1),
    labels = c("Absent", "Present"),
    title = "Domain"
  ),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  rect_gp = gpar(col = "dimgrey", lwd = 0.5)
)


print(ht_57_ps_shared_domains)




####################### Find DDI between interacting partners and intersecting genes
### find Domain-domain interactions between yajL domain and domains of ps interaction proteins
library(tidyverse)
library(data.table)

#files needed for analysis
# hmmer_57_yajl_ps_prntrs_domtblout
# domtblout_list_htp_ps_intersections$yajL_curated_pfam_results


## 1 Load and parse DOMINE database file
# Parse interaction file
get_domine_interactions <- function(file_path) {
  # First detect number of columns
  first_line <- readLines(file_path, n = 1)
  n_cols <- length(strsplit(first_line, "\\|")[[1]])
  
  # Define all possible column names
  all_possible_cols <- c("domain1", "domain2", "PDB", "3did", "Pfam", "SCOP", 
                         "DIP", "EC", "Kegg", "Reactome", "HPRD", "BIND", "IntAct",
                         "MINT", "UniProt", "Experiment", "Prediction", "HC")
  
  # Read with correct number of columns
  interactions <- fread(file_path, sep = "|", header = FALSE) %>%
    setNames(all_possible_cols[1:n_cols]) %>%
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
  
  # Filter for reliable evidence (adjust as needed)
  interactions <- interactions %>%
    filter(PDB == 1 | `3did` == 1 | Pfam == 1 | Experiment == 1 | HC %in% c("HC", "MC", "LC"))
  
  return(interactions)
}


# Load DOMINE database
domine_db_interactions <- get_domine_interactions("C:\\Users\\User\\Desktop\\Protein interactions\\Data collection\\Get common proteins - positive, negative, datasets\\domine-tables-2.0\\INTERACTION.txt")


## 2 Identify DDIs between your gene and partners

# Function to parse domain descriptions from HMMER results
get_domain_info <- function(domtblout) {
  domtblout %>%
    dplyr::select(domain_id = target_accession, 
                  domain_name = target_name) %>%
    mutate(domain_id = str_remove(domain_id, "\\..+$")) %>%  # remove the digits after the dot, eg "PF01965.30" → "PF01965"
    distinct()
}


# Function to analyse DDIs
find_domain_interactions <- function(gene_domains, partner_domains, domine_interactions) {
  
  # Process gene domains
  gene_domain_info <- get_domain_info(gene_domains)
  
  # Process partner domains with multiple ID formats
  # Process partner data - create all domain variants
  partner_data <- partner_domains %>%
    dplyr::select(
      partner_protein_name = query_name,
      partner_domain_id = target_accession,
      partner_domain_name = target_name
    ) %>%
    mutate(partner_domain_id = str_remove(partner_domain_id, "\\..+$")) %>%  # Clean domain ID
    dplyr::select(partner_protein_name, 
                  partner_domain_id, 
                  partner_domain_name) %>%
    distinct()
  
  # Create all possible domain pairs
  domain_pairs <- crossing(
    gene_domain_id = gene_domain_info$domain_id,
    partner_protein_name = unique(partner_data$partner_protein_name)
  ) %>%
    left_join(partner_data, by = "partner_protein_name", relationship = "many-to-many")
  
  # Find interactions with all evidence types
  interactions <- domain_pairs %>%
    left_join(
      domine_interactions,
      by = c("gene_domain_id" = "domain1", "partner_domain_id" = "domain2"),
      relationship = "many-to-many"
    ) %>%
    left_join(
      domine_interactions,
      by = c("gene_domain_id" = "domain2", "partner_domain_id" = "domain1"),
      suffix = c(".a", ".b")
    ) %>%
    # Classify the evidence and its strength (weak, high)
    mutate(
      # Interaction exists if ANY evidence is non-NA/non-zero
      interaction = ifelse(
        (!is.na(PDB.a) & PDB.a == 1) | (!is.na(PDB.b) & PDB.b == 1) |
          (!is.na(`3did.a`) & `3did.a` == 1) | (!is.na(`3did.b`) & `3did.b` == 1) |
          (!is.na(Experiment.a) & Experiment.a == 1) | (!is.na(Experiment.b) & Experiment.b == 1) |
          HC.a %in% c("HC","MC","LC") | HC.b %in% c("HC","MC","LC") |
          (!is.na(Pfam.a) & Pfam.a == 1) | (!is.na(Pfam.b) & Pfam.b == 1),
        TRUE, FALSE
      ),
      
      # Evidence classification (priority order)
      evidence = case_when(
        (!is.na(PDB.a) & PDB.a == 1) | (!is.na(PDB.b) & PDB.b == 1) ~ "High (PDB)",
        (!is.na(`3did.a`) & `3did.a` == 1) | (!is.na(`3did.b`) & `3did.b` == 1) ~ "High (3DID)",
        (!is.na(Experiment.a) & Experiment.a == 1) | (!is.na(Experiment.b) & Experiment.b == 1) ~ "High (Experimental)",
        HC.a == "HC" | HC.b == "HC" ~ "High (HC prediction)",
        HC.a == "MC" | HC.b == "MC" ~ "Medium (MC prediction)",
        HC.a == "LC" | HC.b == "LC" ~ "Low (LC prediction)",
        (!is.na(Pfam.a) & Pfam.a == 1) | (!is.na(Pfam.b) & Pfam.b == 1) ~ "Low (Pfam only)",
        TRUE ~ "No evidence"
      )
    ) %>%
    filter(interaction) %>%
    left_join(gene_domain_info, by = c("gene_domain_id" = "domain_id")) %>%
    dplyr::select(
      gene_domain_id,
      gene_domain_name = domain_name,  # CORRECTED: Use the right column name
      partner_protein_name,
      partner_domain_id,
      partner_domain_name,
      evidence
    ) %>%
    distinct()
  
  return(interactions)
}


# Find interactions
ddi_yajl_ps_results <- find_domain_interactions(
  gene_domains = domtblout_list_htp_ps_intersections$yajL_curated_pfam_results,
  partner_domains = hmmer_57_yajl_ps_prntrs_domtblout,
  domine_interactions = domine_db_interactions
)

# View results
print(ddi_yajl_ps_results)

# #save the results
# library(writexl)
# 
# write_xlsx(
#   ddi_yajl_ps_results,  # Your results tibble
#   path = "DDI_yajl_ps_results.xlsx"  # Output filename
# )





####################### Find Domain Families across interacting partners to identify patterns
#files needed for analysis
# hmmer_57_yajl_ps_prntrs_domtblout


## 1 Load and parse Pfam metadata file
# Parse domain info file
library(stringr)

parse_pfam_dat <- function(file) {
  data <- readLines(file)
  entries <- split(data, cumsum(data == "//"))
  
  map_dfr(entries, ~{
    # Extract all GF fields
    gf_lines <- str_subset(.x, "^#=GF")
    fields <- str_match(gf_lines, "^#=GF\\s+(\\w+)\\s+(.*)")[, 2:3]
    
    # Handle duplicate field names by collapsing values
    field_list <- split(fields[,2], fields[,1])
    field_list <- map(field_list, ~paste(.x, collapse = "|"))
    
    # Convert to tibble with automatic name repair
    as_tibble(field_list, .name_repair = "unique")
  })
}

# Load pfam file
pfam_metadata <- parse_pfam_dat("C:\\Users\\User\\Desktop\\Protein interactions\\Data collection\\Get common proteins - positive, negative, datasets\\Pfam-A.hmm.dat\\Pfam-A.hmm.dat") %>%
  dplyr::select(ID, AC, DE, TP, CL) %>%
  distinct(AC, .keep_all = TRUE) %>%  # Remove duplicate entries
  mutate(
    Pfam_ID = str_extract(AC, "PF\\d+"),
    Type = case_when(
      TP == "Domain" ~ "True domain",
      TP == "Family" ~ "Protein family",
      TP == "Repeat" ~ "Structural repeat",
      TRUE ~ "Other"
    )
  )


## Extract Pfam accessions from the domain table
# Extract just the PFxxxxx part from target_accession
hmmer_57_yajl_ps_prntrs_domtblout <- hmmer_57_yajl_ps_prntrs_domtblout %>%
  mutate(Pfam_ID = str_extract(target_accession, "PF\\d+"))


# Merge with Pfam metadata to get CL (family) information
domains_with_family_57_ps_prntrs <- hmmer_57_yajl_ps_prntrs_domtblout %>%
  left_join(pfam_metadata %>% dplyr::select(Pfam_ID, CL), by = "Pfam_ID")


# Filter missing CL values for domains that have no CL assigned to 
domains_filtered_57_ps_prntrs <- domains_with_family_57_ps_prntrs %>%
  filter(!is.na(CL))  # Remove rows with missing CL

# Create CL × Pfam binary matrix
cl_domain_matrix_57_ps_prntrs <- domains_filtered_57_ps_prntrs %>%
  distinct(target_accession, CL) %>%    # Ensure unique domain-CL pairs
  mutate(in_cluster = 1) %>%           # Binary indicator
  pivot_wider(
    names_from = CL,                   # CL clusters as columns
    values_from = in_cluster,          # Fill with 1/0
    values_fill = 0                    # Absent = 0
  ) %>%
  column_to_rownames("target_accession") %>%  # Domains as row names
  as.matrix()  



library(ComplexHeatmap)
cl_domain_57_ps_prntrs <- Heatmap(
  cl_domain_matrix_57_ps_prntrs,
  col = c("0" = "white", "1" = "darkblue"),
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",        
  column_names_side = "bottom",
  column_title = "Domain-Family Associations",
  row_names_gp = gpar(fontsize = 6),      # Domain label size
  column_names_gp = gpar(fontsize = 8),   # Protein header size
  column_names_rot = 90,                  # Rotated protein names
  heatmap_legend_param = list(
    at = c(0, 1),
    labels = c("Absent", "Present"),
    title = "Domain"
  ),
  cluster_rows = TRUE,           # No row clustering
  cluster_columns = TRUE,        # No column clustering
  show_row_dend = FALSE,          # Explicitly hide row dendrogram
  show_column_dend = FALSE,       # Explicitly hide column dendrogram
  row_order = 1:nrow(cl_domain_matrix_57_ps_prntrs),    # Maintain original row order
  column_order = 1:ncol(cl_domain_matrix_57_ps_prntrs),
  rect_gp = gpar(col = "grey", lwd = 0.5)
)

print(cl_domain_57_ps_prntrs)

# save the heatmap
png("cl_domain_57_ps_prntrs_heatmap.png",
    width = 1600,         # 300 dpi
    height = 2300,
    res = 200)
draw(cl_domain_57_ps_prntrs, padding = unit(c(2, 2, 10, 2), "mm"))
dev.off()
























