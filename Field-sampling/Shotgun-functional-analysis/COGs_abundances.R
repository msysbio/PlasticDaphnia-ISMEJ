"
Deseq2: Compare COGs composition across samples
"

library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)  # For theme_pubr()
library(vegan)

#Load data and packages:
counts = read.csv('/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/Doctorate/Field/Shotgun/14_All_together_datasets/COGs/COGs_all_samples_raw_counts.tsv',row.names = 1,sep="\t")
meta = read.csv('/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/Doctorate/Field/Shotgun/metadata.csv',row.names = 1, header=TRUE)
#colnames(meta)[2] = "condition"
meta <- meta[colnames(counts),]

##remove empty COGs
counts_tidy <- counts %>%
  dplyr::filter(rowSums(.) > 0) #proportion

#take proportions
counts_tidy = data.frame(t(counts_tidy))
counts_total = rowSums(counts_tidy)
proportions_counts = sweep(counts_tidy, 1, counts_total, "/")

bray <- proportions_counts %>%
  vegdist(method = "bray")

#PERMANOVAS
##Between bacterioplankton and Daphnia

# Perform betadisper to assess homogeneity of variances
dispersion_test <- betadisper(bray, meta$Type)

# To see the results
print(dispersion_test)

# To visualize the dispersion
plot(dispersion_test)
anova(dispersion_test) #0.01427 * #Beta-disper is significant. Permanova cannot be performed.

##Between High vs Low - Daphnia
meta_daphnia = meta %>%
  filter(Type=="Daphnia")

meta <- subset(meta, !meta$Location %in% c("WWT Leuven", "WWTP Olsene", "WWTP Harelbeke")))

proportions_counts_daphnia <- proportions_counts %>%
  filter(rownames(proportions_counts) %in% rownames(meta_daphnia))

bray_daphnia <- proportions_counts_daphnia %>%
  vegdist(method = "bray")

# Perform betadisper to assess homogeneity of variances
dispersion_test <- betadisper(bray_daphnia, meta_daphnia$Category)

# To see the results
print(dispersion_test)

# To visualize the dispersion
plot(dispersion_test)
anova(dispersion_test) #0.0825 .

#adonis
# Performing PERMANOVA using the adonis function from the vegan package
permanova_result <- adonis2(bray_daphnia~Category,data=meta_daphnia, method = "bray")

# Viewing the results
print(permanova_result) #0.003 **

# Constrained PCoA (db-RDA) using the 'Category' variable
PCoA_constrained <- capscale(bray_daphnia ~ Category+Location, data = meta_daphnia)

pcoa_scores <- scores(PCoA_constrained, display = "sites")
pcoa_df <- as.data.frame(pcoa_scores)
Daphnia_meta= meta_daphnia[rownames(pcoa_df), ]
Daphnia_meta$Location[Daphnia_meta$Location == "ST. Donatus"] <- "St. Donatus"
Daphnia_meta$Location[Daphnia_meta$Location == "Bourgoyen "] <- "Bourgoyen"
pcoa_df$Category = Daphnia_meta$Category
pcoa_df$Location = Daphnia_meta$Location

CCA1_variance=round(PCoA_constrained$CCA$eig[1]/sum(PCoA_constrained$CCA$eig),2)*100
CCA2_variance=round(PCoA_constrained$CCA$eig[2]/sum(PCoA_constrained$CCA$eig),2)*100

# Define shapes for each location
location_shapes <- c("LRV" = 16, "St. Donatus" = 17, "Kluizen" = 18, "Evangelie" = 15, "Blauwe poort" = 19,"Bourgoyen"=7,
                     "Rotselaar"=8,"Citadell"=10, "MVR"=5)
pcoa_df$Shape <- location_shapes[pcoa_df$Location]  # Assign shapes based on Location

# Define colors for Category
category_colours <- c("Low" = "blue", "High" = "red")
pcoa_df$Colour <- category_colours[pcoa_df$Category]

# Create an empty base plot for the constrained ordination
plot(pcoa_df$CAP1, pcoa_df$CAP2, type = "n", main = expression(italic("Daphnia") ~ "microbiome"),
     xlab = paste("CA1: ", CCA1_variance, ".0%",sep=""),
     ylab = paste("CA2: ", CCA2_variance, ".0%",sep=""))

# Add points with color for Category and shape for Location
points(pcoa_df$CAP1, pcoa_df$CAP2, col = pcoa_df$Colour, pch = pcoa_df$Shape, cex = 2.5)

# Draw convex hulls for each Category
unique_categories <- unique(pcoa_df$Category)
for (category in unique_categories) {
  # Subset data for each category
  cat_points <- pcoa_df[pcoa_df$Category == category, ]
  
  # Calculate the convex hull for the points
  hull_indices <- chull(cat_points$CAP1, cat_points$CAP2)
  hull_indices <- c(hull_indices, hull_indices[1])  # Close the hull by repeating the first point
  
  # Draw the convex hull
  lines(cat_points$CAP1[hull_indices], cat_points$CAP2[hull_indices], col = category_colours[category], lwd = 2)
}

# Add legends for clarity
legend("topright", legend = names(category_colours), col = category_colours, pch = 16, title = "", cex = 0.8, bty = "n")


##Between High vs Low - Bacterioplankton
meta_bac = meta %>%
  filter(Type=="Bacterioplankton")

meta_bac= subset(meta_bac, !meta_bac$Location %in% c("WWTP Leuven", "WWTP Olsene", "WWTP Harelbeke","WWT Leuven"))

proportions_counts_bac <- proportions_counts %>%
  filter(rownames(proportions_counts) %in% rownames(meta_bac))

bray_bac <- proportions_counts_bac %>%
  vegdist(method = "bray")

# Perform betadisper to assess homogeneity of variances
dispersion_test <- betadisper(bray_bac, meta_bac$Category)

# To see the results
print(dispersion_test)

# To visualize the dispersion
plot(dispersion_test)
anova(dispersion_test) #0.01433 *

#Make a constrained PCoA plot
# Constrained PCoA (db-RDA) using the 'Category' variable + Location
PCoA_constrained <- capscale(bray_bac ~ Category+Location, data = meta_bac)

pcoa_scores <- scores(PCoA_constrained, display = "sites")
pcoa_df <- as.data.frame(pcoa_scores)

pcoa_df$Category = meta_bac$Category
pcoa_df$Location = meta_bac$Location

CCA1_variance=round(PCoA_constrained$CCA$eig[1]/sum(PCoA_constrained$CCA$eig),2)*100
CCA2_variance=round(PCoA_constrained$CCA$eig[2]/sum(PCoA_constrained$CCA$eig),2)*100

pcoa_df$Location[pcoa_df$Location == "De Gavers"] <- "De gavers"

# Define shapes for each location
location_shapes <- c("LRV" = 16, "St. Donatus" = 17, "Kluizen" = 18, "Evangelie" = 15, "Blauwe poort" = 19,"Bourgoyen"=7,
                     "Rotselaar"=8,"Citadell"=10, "MVR"=5,"De gavers"=6,"Donk"=7,"Kulak"=9)
pcoa_df$Shape <- location_shapes[pcoa_df$Location]  # Assign shapes based on Location

# Define colors for Category
category_colours <- c("Low" = "blue", "High" = "red")
pcoa_df$Colour <- category_colours[pcoa_df$Category]

# Create the base plot without axes and without the main title
plot(pcoa_df$CAP1, pcoa_df$CAP2, type = "n", main = "",
     xlab = paste("CA1: ", round(CCA1_variance, 1), ".0%",sep=""),
     ylab = paste("CA2: ", round(CCA2_variance, 1), ".0%",sep=""),
     xlim = c(-1, 1.1))  # Extend xlim slightly to ensure "1" is visible

# Add a non-bold title using mtext
mtext("Bacterioplankton", side = 3, line = 1, cex = 1.2, font = 1)

# Add points with color for Category and shape for Location
points(pcoa_df$CAP1, pcoa_df$CAP2, col = pcoa_df$Colour, pch = pcoa_df$Shape, cex = 2.5)

# Draw convex hulls for each Category
unique_categories <- unique(pcoa_df$Category)
for (category in unique_categories) {
  # Subset data for each category
  cat_points <- pcoa_df[pcoa_df$Category == category, ]
  
  # Calculate the convex hull for the points
  hull_indices <- chull(cat_points$CAP1, cat_points$CAP2)
  hull_indices <- c(hull_indices, hull_indices[1])  # Close the hull by repeating the first point
  
  # Draw the convex hull
  lines(cat_points$CAP1[hull_indices], cat_points$CAP2[hull_indices], col = category_colours[category], lwd = 2)
}

# Add legends for clarity
legend("topright", legend = names(category_colours), col = category_colours, pch = 16, title = "", cex = 0.8, bty = "n", y.intersp = 0.8)

###End of permanovas

#Divide into Bacterioplankton / Daphnia

###DAPHNIA
ps=phyloseq(otu_table(counts,taxa_are_rows = TRUE), sample_data(meta))
ps_dap=subset_samples(ps,Type == "Daphnia")
counts_dap=data.frame(otu_table(ps_dap))
meta_dap=data.frame(sample_data(ps_dap))

dds <- DESeqDataSetFromMatrix(countData = counts_dap,
                              colData = meta_dap,
                              design = ~ Category)
dds <- DESeq(dds)
res <- results(dds)
res

#Remove low count COGs
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
res = results(dds)
res05 = results(dds,alpha=0.05) #p-value <0.05
summary(res05)

"
out of 16960 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 6093, 36%
LFC < 0 (down)     : 2943, 17%
outliers [1]       : 0, 0%
low counts [2]     : 699, 4.1%
(mean count < 1)

"

###BACTEROLANKTON
ps=phyloseq(otu_table(counts,taxa_are_rows = TRUE), sample_data(meta))
ps_bac=subset_samples(ps,Type == "Bacterioplankton")
counts_bac=data.frame(otu_table(ps_bac))
meta_bac=data.frame(sample_data(ps_bac))

dds <- DESeqDataSetFromMatrix(countData = counts_bac,
                              colData = meta_bac,
                              design = ~ Category)
dds <- DESeq(dds)
res <- results(dds)
res

#Remove low count COGs
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
res = results(dds)
res05 = results(dds,alpha=0.05) #p-value <0.05
summary(res05)

"
out of 54642 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 5730, 10%
LFC < 0 (down)     : 7548, 14%
outliers [1]       : 0, 0%
low counts [2]     : 7554, 14%
(mean count < 2)

"