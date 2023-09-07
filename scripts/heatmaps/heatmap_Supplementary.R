# making a heatmap for read depth, coverage and SNPs to showcase differences in host read filtering with or without EVEs

library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(tidyr)

# home directory
setwd("C:/Users/nadja/Documents/R/trees_for_quaranja_new/MAPPING")

# Set the path to the folder containing the text files
TR_directory <- "C:/Users/nadja/Documents/R/trees_for_quaranja_new/MAPPING"

# Get the list of text files in the folder
file_list <- list.files(TR_directory, pattern = "*__GE_coverage.txt", full.names = TRUE)

# Read all text files into a list of data frames
data_list <- map(file_list, read_delim, delim = "\t", show_col_types = FALSE)

# Combine all data frames into a single data frame
TR_GE <- bind_rows(data_list)

# Add suffix "_TR" to the headers of TR_df
names(TR_GE)<- paste0(names(TR_GE), "_GE")
TR_GE <- TR_GE %>%
  mutate(ID = sub("__GE$", "", `#rname_GE`))
TR_GE <- TR_GE %>%
  mutate(across(everything(), ~str_replace_all(., "-", "_")))

# Set the path to the folder containing the text files
GN_directory <- "C:/Users/nadja/Documents/R/trees_for_quaranja_new/MAPPING/"

# Get the list of text files in the folder
file_list <- list.files(GN_directory, pattern = "*__GWE_coverage.txt", full.names = TRUE)

# Read all text files into a list of data frames
data_list <- map(file_list, read_delim, delim = "\t", show_col_types = FALSE)

# Combine all data frames into a single data frame
GWE_df <- bind_rows(data_list)

names(GWE_df)<- paste0(names(GWE_df), "_GWE")

GWE_df <- GWE_df %>%
  mutate(ID = sub("__GWE$", "", `#rname_GWE`))
GWE_df <- GWE_df %>%
  mutate(across(everything(), ~str_replace_all(., "-", "_")))

# Set the path to the folder containing the text files
TE_directory <- "C:/Users/nadja/Documents/R/trees_for_quaranja_new/MAPPING"

# Get the list of text files in the folder
file_list <- list.files(TE_directory, pattern = "*__TE_coverage.txt", full.names = TRUE)

# Read all text files into a list of data frames
data_list <- map(file_list, read_delim, delim = "\t", show_col_types = FALSE)

# Combine all data frames into a single data frame
TE_df <- bind_rows(data_list)
names(TE_df)<- paste0(names(TE_df), "_TE")
TE_df <- TE_df %>%
  mutate(ID = sub("__TE$", "", `#rname_TE`))
TE_df <- TE_df %>%
  mutate(across(everything(), ~str_replace_all(., "-", "_")))


# Set the path to the folder containing the text files
TE_directory <- "C:/Users/nadja/Documents/R/trees_for_quaranja_new/MAPPING"

# Get the list of text files in the folder
file_list <- list.files(TE_directory, pattern = "*__TWE_coverage.txt", full.names = TRUE)

# Read all text files into a list of data frames
data_list <- map(file_list, read_delim, delim = "\t", show_col_types = FALSE)

# Combine all data frames into a single data frame
TWE_df <- bind_rows(data_list)
names(TWE_df)<- paste0(names(TWE_df), "_TWE")
TWE_df <- TWE_df %>%
  mutate(ID = sub("__TWE$", "", `#rname_TWE`))
TWE_df <- TWE_df %>%
  mutate(across(everything(), ~str_replace_all(., "-", "_")))

# Merge the dataframes based on "ID"
merged_df <- merge(TR_GE, GWE_df, by = "ID", all = TRUE)
merged_df <- merge(merged_df, TE_df, by = "ID", all = TRUE)
merged_df <- merge(merged_df, TWE_df, by = "ID", all= TRUE)


# read depth

meandepth_df <- select(merged_df, ID, meandepth_GE, meandepth_GWE, meandepth_TE, meandepth_TWE)
meandepth_df[is.na(meandepth_df)] <- 0
meandepth_df <- meandepth_df %>%
  mutate_at(vars(-1), as.numeric)

# normalize to highest value per row
meandepth_df[2:5] <-  meandepth_df[2:5]/do.call(pmax, meandepth_df[2:5])

# Transform the dataframe for geom_tile plot
df_transformed <- meandepth_df %>%
  pivot_longer(cols = starts_with("meandepth"), names_to = "meandepth_type", values_to = "meandepth_value") %>%
  mutate(z=as.numeric(scale(meandepth_value)))

df_transformed[df_transformed == 0] <- NA

#df_transformed$Sqrt.meandepth <- sqrt(df_transformed$meandepth_value)

meandepth_plot <- ggplot(data = df_transformed, mapping = aes(x = meandepth_type,
                                            y = ID, fill = meandepth_value)) +
  geom_tile() + #color = "black") +
  xlab(label = "mean sequencing depth") +
  #ylab(label = "Segment") + xlab(label = "SRA sample") +
  #scale_fill_gradient(low = "#d9d9d9", high = "#cb181d", na.value = 'white', name = "average cov\nlength %" ) +
  # keep them squared
  coord_fixed(0.1) +
  scale_fill_gradient(high = "#deebf7", low = "#012345", name= "normalized\n mean depth") +
  scale_x_discrete(position = "top", labels=c("GE", "GWE", "TE", "TWE")) +
#scale_y_continuous(minor_breaks = c(1.7, 2.3, 4.1)) +
# change the width and the height of the legend color bar
#guides(fill = guide_colourbar(barwidth = 0.5,
#barheight = 5, title = "Coverage %")) +
#theme_bw() 
  theme(axis.text.y = element_text(color = "grey20", size = 7, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        panel.background = element_rect(fill = "white"), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, size = 8, hjust = 0), axis.title.y =element_blank()) 

meandepth_plot


# coverage lengths

table <- read.csv("all_sequences.csv")
Xlength <- lengths(regmatches(table$Sequence, gregexpr("N", table$Sequence)))
totallength <- nchar(table$Sequence)
Percentage <- 100-(Xlength/totallength*100)
table$Percentage <- Percentage

table <- table %>%
  mutate(across(everything(), ~str_replace_all(., "-", "_")))
# Separate the "Name" column into two columns: "Name" and "coverage_type"
table <- separate(table, col = "Name", into = c("Name", "coverage_type"), sep = "__", remove = FALSE)


# Separate the "Name" column into two columns: "Name" and "coverage_type"
#table <- separate(table, col = "Name", into = c("Name", "coverage_type"), sep = "__", remove = FALSE)

df_merged <- table %>%
  pivot_wider(
    names_from = coverage_type,
    values_from = c(Sequence, Percentage),
    names_sep = "_"
  )

df_merged <- select(df_merged, Name, Percentage_GE, Percentage_GWE, Percentage_TE, Percentage_TWE)
df_merged <- df_merged %>%
  mutate_at(vars(2:5), as.numeric)
df_merged[is.na(df_merged)] <- 0
df_merged <- df_merged %>%
  mutate_at(vars(-1), as.numeric)

# normalize to highest value per row
df_merged[2:5] <-  df_merged[2:5]/do.call(pmax, df_merged[2:5])

# Transform the dataframe for geom_tile plot
df_transformed <- df_merged %>%
  pivot_longer(cols = starts_with("Percentage"), names_to = "coverage_type", values_to = "coverage_value") 
#df_transformed[df_transformed == 0] <- NA

coverage_plot <- ggplot(data = df_transformed, mapping = aes(x = coverage_type,
                                                              y = Name, fill = coverage_value)) +
  geom_tile() + #color = "black") +
  xlab(label = "coverage analysis") +
  #ylab(label = "Segment") + xlab(label = "SRA sample") +
  #scale_fill_gradient(low = "#d9d9d9", high = "#cb181d", na.value = 'white', name = "average cov\nlength %" ) +
  # keep them squared
  coord_fixed(0.1) +
  scale_fill_gradient(low = "darkred", high = "#fcbba172", name= "normalized\n coverage length") +
  scale_x_discrete(position = "top",labels=c("GE", "GWE", "TE", "TWE")) +
  #scale_y_continuous(minor_breaks = c(1.7, 2.3, 4.1)) +
  # change the width and the height of the legend color bar
  #guides(fill = guide_colourbar(barwidth = 0.5,
  #barheight = 5, title = "Coverage %")) +
  #theme_bw() 
  theme(axis.text.y = element_text(color = "grey20", size = 7, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        panel.background = element_rect(fill = "white"), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45), axis.title.y =element_blank()) 

coverage_plot

# identity plot 

# Load the required libraries
library(Biostrings)
library(dplyr)
library(tidyr)

# Function to count SNPs between two sequences excluding Ns
count_snps <- function(ref_sequence, query_sequence) {
  ref_bps <- as.character(ref_sequence)
  query_bps <- as.character(query_sequence)
  
  # Exclude Ns from the SNP count
  #ref_bps_no_N <- ref_bps[ref_bps != "N"]
  #query_bps_no_N <- query_bps[query_bps != "N"]
  
  #count <- sum(ref_bps_no_N != query_bps_no_N, na.rm = TRUE)
  #return(count)
  count <- sum(ref_bps != query_bps, na.rm = TRUE)
  return(count)
}

# Function to process each alignment file
process_alignment_file <- function(file_path) {
  # Read the alignment file
  alignment <- readDNAStringSet(file_path, format = "fasta")
  
  # Identify the reference sequence with "__TWE" in the header
  reference_name <- grep("__GE", names(alignment), value = TRUE, fixed = TRUE)
  
  if (length(reference_name) != 1) {
    reference_name <- NA
  }
  
  # Initialize a list to store the SNP counts
  snp_counts <- list()
  
  # Loop through each sequence and count SNPs to the reference
  for (name in names(alignment)) {
    sequence <- as.character(alignment[name])
    if (is.na(reference_name)) {
      snp_count <- NA
    } else {
      snp_count <- count_snps(alignment[reference_name], sequence)
    }
    snp_counts[[name]] <- snp_count
  }
  
  # Combine SNP counts into a data frame
  result_df <- data.frame(Sequence = names(snp_counts), SNPs = unlist(snp_counts))
  
  return(result_df)
}

# List of file paths for alignment files
alignment_files <- list.files(path = "C:/Users/nadja/Documents/R/trees_for_quaranja_new/MAPPING",
                              pattern = "alignment\\.fasta$", full.names = TRUE)

# Process each alignment file and bind the results into a final table
final_table <- lapply(alignment_files, process_alignment_file) %>%
  bind_rows()
#write.table(final_table, file = "table.txt", row.names = FALSE, quote = FALSE)
final_table <- read.table("table.txt", header = TRUE)
# Separate the "Name" column into two columns: "Name" and "coverage_type"
final_table <- separate(final_table, col = "Sequence", into = c("Sequence", "coverage_type"), sep = "__", remove = FALSE)


# Print the final table
final_table$SNPs <- as.numeric(final_table$SNPs)
# Convert SNPs to a factor with three discrete categories
final_table$SNPs <- cut(final_table$SNPs, breaks = c(-Inf, 0.5, 1.5, Inf), labels = c(0, 1, 2))

colors <- c("white", "darkblue", "darkred")

SNP_plot <- ggplot(data = final_table, mapping = aes(x = coverage_type, y = Sequence, fill = SNPs)) +
  geom_tile(color = "black") +
  xlab(label = "SNP analysis") +
  coord_fixed(0.1) +
  scale_fill_manual(values = colors, labels = c("0 SNPs", "1 SNP", "2 SNPs")) +
  scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(color = "grey20", size = 7, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.title.y = element_blank())

print(SNP_plot)


library(patchwork)
meandepth_plot + coverage_plot +  SNP_plot
