"
PlasticDaphnia: Nicola's exposure
Bacterioplankton vs Microbiome
"
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggpubr)
library(rgl)
library(BiocManager)
library(microbiome)
library(knitr)
#install.packages("patchwork")
library(patchwork)

"Bacterioplankton + Microbiome"

#Load phyloseq (non-rarefied)
ps_dap = readRDS("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Exposure-Niki/tables/phyloseq/decontam/ps_daphnia_decontam.rds")
ps_bac = readRDS("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Exposure-Niki/tables/phyloseq/decontam/ps_bacterioplankton_decontam.rds")

#Merge phyloseq
merged_phyloseq <- merge_phyloseq(ps_dap, ps_bac)
#Rarefy
sort(sample_sums(merged_phyloseq), decreasing=FALSE)
# Set your chosen depth
set_depth <- 1563

# Perform rarification
#merged_phyloseq_rarified <- rarefy_even_depth(merged_phyloseq, sample.size = set_depth, rngseed = 1, replace = FALSE)
#saveRDS(merged_phyloseq_rarified,"/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Exposure-Niki/tables/phyloseq/decontam/merged_phyloseq_rarified.rds")
merged_phyloseq_rarified = readRDS("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Exposure-Niki/tables/phyloseq/decontam/merged_phyloseq_rarified.rds")
meta = data.frame(sample_data(merged_phyloseq_rarified))
#write.csv(meta,"/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Exposure-Niki/tables/decontaminated_set/merged_meta.csv")
#meta = read.csv("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Exposure-Niki/tables/decontaminated_set/merged_meta.csv",row.names = 1)
sample_data(merged_phyloseq) = meta

#Post-exposure control (no MPs) daphnia vs bacterioplankton (no MPs)

subset_phyloseq <- subset_samples(merged_phyloseq_rarified, 
                                  (Category == "Daphnia: Post-exposure" & Plastic == "C" & Pond!="H2O") |
                                    (Category == "Bacterioplankton: Post-exposure" & Plastic == "C" & Pond!="H2O"))

meta = data.frame(sample_data(subset_phyloseq))

bray <- vegdist(data.frame(t(otu_table(subset_phyloseq))), method="bray") #creating dissimilarity matrix with vegdist using bray curtis
nmds_result <- metaMDS(bray, trymax = 100,k=2)  # stress = 0.1785537 


#Permanova

# Run PERMANOVA using adonis
permanova_result <- adonis2(bray~ Category, data = meta, method = "bray")
# View the results
print(permanova_result) #0.019 *

#Check assumptions of constant variance
dispersion_test <- betadisper(bray, meta$Category)
plot(dispersion_test)
permutest(dispersion_test,permutations = 999) #0.727

#Richness

tab <-microbiome::alpha(subset_phyloseq, index = "all")
tab$Category = meta$Category


###Chao1
p.chao1 <- boxplot_alpha(subset_phyloseq, 
                         index = "Chao1",
                         x_var = "Category",
                         fill.colors = c("Daphnia: Post-exposure" = "yellow", 
                                         "Bacterioplankton: Post-exposure" = "green"))


p.chao1 <- p.chao1 + theme_pubr() + 
  labs(x="", y="Chao1") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.position = "none") # This line removes the legend

p.chao1
print(t.test(chao1~Category, data=tab)) # p-value = 0.2637


#####Shannon

p.shannon <- boxplot_alpha(subset_phyloseq, 
                           index = "diversity_shannon",
                           x_var = "Category",
                           fill.colors = c("Daphnia: Post-exposure" = "yellow", 
                                           "Bacterioplankton: Post-exposure" = "green"))


p.shannon <- p.shannon + theme_pubr() + 
  labs(x="", y="Shannnon") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.position = "none") # This line removes the legend

p.shannon
print(t.test(diversity_shannon~Category, data=tab)) # p-value = 0.2106


#####Pielou

p.Pielou <- boxplot_alpha(subset_phyloseq, 
                          index = "evenness_pielou",
                          x_var = "Category",
                          fill.colors = c("Daphnia: Post-exposure" = "yellow", 
                                          "Bacterioplankton: Post-exposure" = "green"))


p.Pielou <- p.Pielou + theme_pubr() + 
  labs(x="", y="Pielou") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.position = "none") # This line removes the legend

p.Pielou
print(t.test(evenness_pielou~Category, data=tab)) # p-value = 0.282

#####Combine in one plot
p_combined <- p.chao1 / p.shannon / p.Pielou
p_combined
# Optionally, if you want to add a common title or adjust spacing, you can do:
p_combined + 
  plot_annotation(title = "Diversity Indices Comparison for Bacterioplantkon samples") +
  plot_layout(guides = 'collect') # This collects and positions the legends together if applicable


#10 TOP TAXA ACROSS HIGH AND LOW-Bacterioplankton

ps_filtered_bac.genus <- tax_glom(ps_filtered_bac, taxrank = "Genus")
bp_NP <- subset_samples(ps_filtered_bac.genus, Category=="Natural pond")
# Top N taxa
N <- 10
top <- names(sort(taxa_sums(bp_NP), decreasing = TRUE))[1:N]

# Calculate relative abundance
ps_filtered_bac.genus.prop <- transform_sample_counts(bp_NP, function(x) x / sum(x) )

# Subset object to top N taxa
ps_filtered_bac.genus.prop.top <- prune_taxa(top, ps_filtered_bac.genus.prop)

ps_filtered_bac.genus.prop.top_df = psmelt(ps_filtered_bac.genus.prop.top) 

ggplot(ps_filtered_bac.genus.prop.top_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Sample", y = "Relative Abundance", title = "Top 10 Taxa Relative Abundance for Natural pond in bacterioplankton") +
  scale_fill_viridis_d(begin = 0.5, direction = -1) 

"End of code"


