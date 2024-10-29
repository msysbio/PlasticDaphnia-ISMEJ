"
PlasticDaphnia: Nicola's exposure
Removal of contaminants
Using external package: microDecon : https://github.com/donaldtmcknight/microDecon
"

#install.packages('/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Exposure-Niki/scripts/microDecon-master', repos = NULL, type = "source")
library(microDecon)
library(phyloseq)
library(vegan)
# Load your data
# Replace these with your actual file paths
otus = read.csv("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Exposure-Niki/tables/SILVA/otus_clean.csv", row.names=1)
taxa = read.csv("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Exposure-Niki/tables/SILVA/tax_clean.csv", row.names=1)
meta = read.csv("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Exposure-Niki/tables/meta.csv", row.names=1)

otus_sorted <- otus[rownames(meta), ]
otus_sorted = data.frame(t(otus_sorted))
ps <- phyloseq(otu_table(otus_sorted, taxa_are_rows = TRUE),
               tax_table(as.matrix(taxa)), sample_data(meta))

##Daphnia
ps_dap = subset_samples(ps, Type == "Daphnia")
otus_dap = data.frame(otu_table(ps_dap))
taxa_dap = data.frame(tax_table(ps_dap))
meta_dap = data.frame(sample_data(ps_dap))

all(rownames(taxa_dap) == rownames(otus_dap))
create_taxa_string <- function(row) {
  labels <- c("K_", "P_", "C_", "O_", "F_", "G_", "S_")
  # Select only the first 7 elements of the row (or however many you need)
  row_selected <- row[1:7]
  # Combine labels with taxa
  taxa_string <- paste0(labels, row_selected, collapse = "; ")
  # Replace '; NA' with empty string in case of NA values
  taxa_string <- gsub("; NA", "", taxa_string)
  # Trim trailing semicolons if any
  taxa_string <- gsub("; $", "", taxa_string)
  return(taxa_string)
}

taxa_naming <- apply(taxa_dap, 1, create_taxa_string)
otus_dap$Taxa = taxa_naming
otus_dap$OTU_ID = rownames(otus_dap)
otus_dap <- otus_dap[c("OTU_ID", setdiff(names(otus_dap), "OTU_ID"))]
colnames(otus_dap)[colnames(otus_dap) == "NC_InsectKit"] = "Blank1"

#Put Blank1 as a second column:
if("Blank1" %in% colnames(otus_dap)) {
  # Move 'Blank1' to the second position
  otus_dap <- otus_dap[c(colnames(otus_dap)[1], "Blank1", setdiff(colnames(otus_dap), c(colnames(otus_dap)[1], "Blank1")))]
} else {
  message("Column 'Blank1' not found in otus_bac")
}

decontaminated <- decon(data = otus_dap,numb.ind=c(23,9,41,53,2),numb.blanks=1,taxa=T,runs=1)
otus_dap_decontam = decontaminated$decon.table
rownames(otus_dap_decontam) =otus_dap_decontam$OTU_ID
otus_dap_decontam = otus_dap_decontam[,-c(1:2)]
otus_dap_decontam = otus_dap_decontam[,-length(colnames(otus_dap_decontam))]

ps_dap_decontam = phyloseq(otu_table(otus_dap_decontam, taxa_are_rows = TRUE),
                           tax_table(as.matrix(taxa_dap)), sample_data(meta_dap))

#Remove chlorophyll
ps_dap_decontam <- subset_taxa(ps_dap_decontam, Order != "Chloroplast")
# Remove mitochondria
ps_dap_decontam <- subset_taxa(ps_dap_decontam, Family != "Mitochondria")

saveRDS(ps_dap_decontam,"/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Exposure-Niki/tables/phyloseq/decontam/ps_daphnia_decontam.rds")
ps_dap_decontam=readRDS("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Exposure-Niki/tables/phyloseq/decontam/ps_daphnia_decontam.rds")
#Sample's depth before deconamination:
rarecurve(t(data.frame(otu_table(ps_dap))), step = 100, col = "blue", cex = 0.7)

#Sample's depth after deconamination:
rarecurve(t(data.frame(otu_table(ps_dap))), step = 100, col = "blue", cex = 0.7)
sort(sample_sums(ps_dap_decontam), decreasing=FALSE)

# Set your chosen depth
set_depth <- 2494  # lowest depth sample

# Perform rarification
ps_dap_rarified <- rarefy_even_depth(ps_dap_decontam, sample.size = set_depth, rngseed = 1, replace = FALSE)
sort(sample_sums(ps_dap_rarified), decreasing=FALSE)

saveRDS(ps_dap_rarified,"/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Exposure-Niki/tables/phyloseq/decontam/ps_rare_daphnia_decontam.rds")

##Bacterioplankton
ps_bac = subset_samples(ps, Type == "Bacterioplankton")
otus_bac = data.frame(otu_table(ps_bac))
taxa_bac = data.frame(tax_table(ps_bac))
meta_bac = data.frame(sample_data(ps_bac))

all(rownames(taxa_bac) == rownames(otus_bac))
taxa_naming <- apply(taxa_bac, 1, create_taxa_string)
otus_bac$Taxa = taxa_naming
otus_bac$OTU_ID = rownames(otus_bac)
otus_bac <- otus_bac[c("OTU_ID", setdiff(names(otus_bac), "OTU_ID"))]
colnames(otus_bac)[colnames(otus_bac) == "NA285"] = "Blank1"

#Put Blank1 as a second column:
if("Blank1" %in% colnames(otus_bac)) {
  # Move 'Blank1' to the second position
  otus_bac <- otus_bac[c(colnames(otus_bac)[1], "Blank1", setdiff(colnames(otus_bac), c(colnames(otus_bac)[1], "Blank1")))]
} else {
  message("Column 'Blank1' not found in otus_bac")
}

decontaminated <- decon(data = otus_bac,numb.ind=c(95,6),numb.blanks=1,taxa=T,runs=1)
otus_bac_decontam = decontaminated$decon.table
rownames(otus_bac_decontam) =otus_bac_decontam$OTU_ID
otus_bac_decontam = otus_bac_decontam[,-c(1:2)]
otus_bac_decontam = otus_bac_decontam[,-length(colnames(otus_bac_decontam))]

#Sample's depth before deconamination:
rarecurve(data.frame(t(otu_table(ps_bac))), step = 100, col = "blue", cex = 0.7)

#Sample's depth after deconamination:
rarecurve(t(otus_bac_decontam), step = 100, col = "blue", cex = 0.7)

ps_bac_decontam = phyloseq(otu_table(otus_bac_decontam, taxa_are_rows = TRUE),
                           tax_table(as.matrix(taxa_bac)), sample_data(meta_bac))

sort(sample_sums(ps_bac_decontam), decreasing=FALSE)

saveRDS(ps_bac_decontam,"/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Exposure-Niki/tables/phyloseq/decontam/ps_bacterioplank_decontam.rds")
ps_bac_decontam=readRDS("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Exposure-Niki/tables/phyloseq/decontam/ps_bacterioplank_decontam.rds")

#Remove chlorophyll
ps_bac_decontam <- subset_taxa(ps_bac_decontam, Order != "Chloroplast")
# Remove mitochondria
ps_bac_decontam <- subset_taxa(ps_bac_decontam, Family != "Mitochondria")


sort(sample_sums(ps_bac_decontam), decreasing=FALSE)
# Set your chosen depth
set_depth <- 1563  # for example, choose based on your sample_depths output

# Perform rarification
ps_bac_rarified <- rarefy_even_depth(ps_bac_decontam, sample.size = set_depth, rngseed = 1, replace = FALSE)

saveRDS(ps_bac_rarified,"/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Exposure-Niki/tables/phyloseq/decontam/ps_rare_bacteroplankton_decontam.rds")

#End of code

