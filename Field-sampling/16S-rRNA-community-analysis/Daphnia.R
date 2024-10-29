"
Daphnia-microbiome: analysis
"

library(phyloseq)
library(vegan)
library(ggplot2)
library(ggpubr)
library(plotly)
library(ggforce)

ps_raw = readRDS("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Miseq/output/dada2/new/ps_raw")
meta <- read.csv("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Miseq/metadata/metadata-basic.csv", row.names=1)

sample_data(ps_raw) = meta

#get daphnia
dap=subset_samples(ps_raw, sample_data(ps_raw)$Type == "Daphnia")
#Remove samples that are not daphnia
samples_to_remove=c("N62","N64")
dap <- subset_samples(dap, !(sample_names(dap) %in% samples_to_remove))
otus = data.frame(otu_table(dap))

# Create a rarefaction curve
rarecurve(otus, step = 100, col = "blue", cex = 0.7)
sort(sample_sums(dap), decreasing=FALSE)

# Set your chosen depth
set_depth <- 8949  # for example, choose based on your sample_depths output
# Perform rarification
dap_rarified <- rarefy_even_depth(dap, sample.size = set_depth, rngseed = 1, replace = FALSE)
saveRDS(dap_rarified, "/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Miseq/output/dada2/new/dap_rarified")
dap_rarified = readRDS("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Miseq/output/dada2/new/dap_rarified")

bray <- vegdist(data.frame(otu_table(dap_rarified)), method="bray") #creating dissimilarity matrix with vegdist using bray curtis
nmds_result <- metaMDS(bray, trymax = 1000, autotransform = FALSE)  # stress = 0.1768015 #GOOD

# Custom color scheme for categories
colors <- c("Daphnia" = "gold", 
            "Artifical pond" = "navy")

# Assuming 'Location' is a column in your metadata
p = NULL
p = plot_ordination(dap_rarified, nmds_result, color = "Category") +
  geom_point(aes(shape = Location), size = 5) +  # Increase point size here
  scale_color_manual(values = colors) +  # Apply custom color scheme
  scale_shape_manual(values = c(0, 1, 2,5,6,15,16,17,18)) +  # Custom shapes for each location
  labs(title = "", 
       x = "NMDS1", 
       y = "NMDS2", 
       color = "Sample Type", 
       shape = "Location") +  # Update labels
  theme_pubclean()+  # Minimal theme
  theme(legend.text = element_text(size = 6),  # Smaller legend text
        legend.title = element_text(size = 5),  # Smaller legend title
        legend.key.size = unit(0.5, "cm"))  # Smaller legend keys

p
#Envfit : add env.paraemeters to the plot

env_meta <- read.csv("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Miseq/metadata/RDA/daphnia/env_parameters.csv", row.names = 1)
# Assuming your OTU table is in a phyloseq object named 'ps'
otu_samples <- sample_names(dap_rarified)
# Assuming your environmental data is in a dataframe named 'env_data'
env_samples <- rownames(env_meta)
# Find common samples
common_samples <- intersect(otu_samples, env_samples)
# Subset environmental data to include only common samples
env_data_subset <- env_meta[common_samples, ]
# Optional: Order the environmental data to match the order in OTU table
env_data_ordered <- env_data_subset[match(otu_samples, rownames(env_data_subset)), ]

#remove gaps
# Simple mean imputation for each column
env_data_filled <- env_data_ordered
for(i in seq_along(env_data_filled)) {
  env_data_filled[[i]][is.na(env_data_filled[[i]])] <- mean(env_data_filled[[i]], na.rm = TRUE)
}

env_data_filled = data.frame(env_data_filled[,-c(1:2)])
#save imputed dataset
#write.csv(env_data_filled, "/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Miseq/metadata/RDA/daphnia/env_parameters_imputed.csv")
env_data_filled = read.csv("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Miseq/metadata/RDA/daphnia/env_parameters_imputed.csv", row.names=1)

# Running envfit
envfit_result <- envfit(nmds_result, env_data_filled, permutations = 999,scaling = "sites")
significant_scores <- data.frame(envfit_result$vectors$arrows[envfit_result$vectors$pvals < 0.05, ])
significant_scores$R = envfit_result$vectors$r[rownames(significant_scores)]

# Assuming 'Location' is a column in your metadata
p = NULL
p = plot_ordination(dap_rarified, nmds_result, color = "Category") +
  geom_point(aes(shape = Location), size = 5) +  # Increase point size here
  scale_color_manual(values = colors) +  # Apply custom color scheme
  scale_shape_manual(values = c(0, 1, 2,5,6,15,16,17,18)) +  # Custom shapes for each location
  labs(title = "", 
       x = "NMDS1", 
       y = "NMDS2", 
       color = "Sample Type", 
       shape = "Location") +  # Update labels
  theme_pubclean()+  # Minimal theme
  theme(legend.text = element_text(size = 6),  # Smaller legend text
        legend.title = element_text(size = 5),  # Smaller legend title
        legend.key.size = unit(0.5, "cm"))  # Smaller legend keys

scaling_constant = 5
# Overlay significant vectors on the base plot
p <- p +
  geom_segment(data = significant_scores, 
               aes(x = 0, y = 0, xend = NMDS1 * R * scaling_constant , yend = NMDS2 * R * scaling_constant), 
               arrow = arrow(type = "closed", length = unit(0.1, "inches"), ends = "last", angle = 20),
               color = "red", alpha = 0.3)  # Adjust transparency
# Use annotate() to add labels to the arrows
p <- p +
  annotate(geom = "text", x = significant_scores$NMDS1 * significant_scores$R * scaling_constant, y = significant_scores$NMDS2 * significant_scores$R * scaling_constant, 
           label = rownames(significant_scores), vjust = -0.5)
p

#PERMANOVA

#Statistics:
##Permanova:
permanova = adonis2(bray ~ data.frame(sample_data(dap_rarified))$Category) #0.001 ***

##Permanova assumption tests:
# test homogeneity of variance
homog_test_cat <- betadisper(bray, data.frame(sample_data(dap_rarified))$Category,bias.adjust=TRUE)
# anova on centroids
disp_category=anova(homog_test_cat) #0.05235 .

#Richness
library(BiocManager)
#BiocManager::install("microbiome")
library(microbiome)
library(knitr)
library(ggpubr)

tab <-microbiome::alpha(dap_rarified, index = "all")
meta <- data.frame(sample_data(dap_rarified)) #Accessing my sample information from the ps object containing rarefied data
tab$Category = meta$Category

p.chao1 <- boxplot_alpha(bap_rarified, 
                           index = "Shannon",
                           x_var = "Category",
                           fill.colors = c("Natural pond" = "gold", 
                                           "Artifical pond" = "navy"))


p.chao1 <- p.chao1 + theme_minimal() + 
  labs(x="", y="Shannon") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))

p.chao1
print(t.test(chao1~Category, data=tab)) # p-value = 0.0000000000002324



