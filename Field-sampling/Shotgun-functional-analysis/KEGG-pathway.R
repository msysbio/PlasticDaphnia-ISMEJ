library(tidyverse)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(pheatmap)
library(grDevices)
library(stringr)  # For text folding

# Load the data
df <- read.csv("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/3_Daphnia_microbiome/KEGG-KOs/KO-pathways-abundances.csv")
df = df %>%
  dplyr::filter(Polymer != "PS")

df$Order_for_heatmap <- as.numeric(as.character(df$Order_for_heatmap))

# Arrange the data frame by 'Order_for_heatmap'
df_ordered <- df %>%
  arrange(Order_for_heatmap) %>%
  mutate(Enzyme = factor(Enzyme, levels = unique(Enzyme)))

# Add KO number to Enzyme and fold descriptions (shortening the gene names)
df_ordered$Enzyme <- paste(df_ordered$KO.number, ":", df_ordered$Enzyme, sep="")

# Shorten or fold the Enzyme descriptions by adding a line break at 30 characters
df_ordered$Enzyme <- str_wrap(df_ordered$Enzyme, width = 30)

# Store enzyme names and pathway information
enzyme_names <- df_ordered$Enzyme
pathway <- df_ordered$Pathway

# Convert data to matrix
data_matrix <- as.matrix(cbind(df_ordered$Artifical.pond, df_ordered$Natural.pond))
data_matrix <- sqrt(data_matrix)
colnames(data_matrix) <- c("High MPs", "Low MPs")
rownames(data_matrix) <- enzyme_names

# Create a color vector for the conditions
pathway_conditions <- df_ordered %>% 
  select(Pathway)

rownames(pathway_conditions) <- df_ordered$Enzyme
annotation_colors <- list(Pathway = c("Styrene degradation" = "darkgreen", 
                                      "Polycyclic aromatic hydrocarbon degradation" = "yellow"))

# Generate heatmap
pheatmap(data_matrix, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         color = colorRampPalette(c("white", "darkred"))(100), 
         display_numbers = TRUE, 
         fontsize_row = 12,  # Adjust row font size
         fontsize_col = 14,  # Increase font size for x-axis labels
         fontsize_number = 12,  # Increase font size inside tiles
         number_color = "black",
         border_color = "black",
         main = "PET degradation pathway",  # Title
         angle_col = 0,  # Make x-axis labels horizontal
         fontface_row = "italic")  # Italicize y-axis labels



