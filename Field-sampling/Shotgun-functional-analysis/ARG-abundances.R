library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

###ARGs abundances comparison
# Load the data
df <- read.csv("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/Doctorate/Publication/tables/Fields shotgun/ARG_abundances.csv",row.names=1)

# Transform data from wide to long format
long_data <- gather(df, Key, Value, -Category)

# Calculate average values for each enzyme family by pond type
average_values <- long_data %>%
  group_by(Category) %>%
  summarize(Average = mean(Value))

# Update the color vector with new labels
colors <- c("Low MPs" = "#6885d0", "High MPs" = "#cb5658")

# Ensure the Category column in long_data has updated labels
long_data$Category <- recode(long_data$Category, "Low" = "Low MPs", "High" = "High MPs")

# Create the boxplot with dot overlay
boxplot_chart <- ggplot(long_data, aes(x=Category, y=Value, fill=Category)) +
  # Boxplot for distribution
  geom_boxplot(outlier.shape = NA, alpha=0.5) +  # Hide outliers in the boxplot itself
  # Jitter to overlay individual points
  geom_jitter(aes(color=Category), width=0.2, size=2, alpha=0.7) +
  # Labels and styling
  labs(y="Average CPM per pond ategory", x="Pond category") +
  scale_fill_manual(values=colors) +  # Fill color for the boxes
  scale_color_manual(values=colors) +  # Color for the dots
  theme_pubr() +  # Optionally, use theme_minimal() for a clean look
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(boxplot_chart)

#AR-genes difference between High and Low MPs ponds
wilcox.test(Average.ARG.abundance ~ Category,data=df,exact = TRUE) #pval=0.02857
"""
Wilcoxon rank sum exact test

data:  Average.ARG.abundance by Category
W = 16, p-value = 0.02857
alternative hypothesis: true location shift is not equal to 0
"""
#Calculate the effect size:
install.packages("rstatix")

# Load necessary library
library(rstatix)

# Effect size
effect_size_result <- wilcox_effsize(Average.ARG.abundance ~ Category, data = df, exact = TRUE)

# Print results
print(effect_size_result)

# Replace 'path_to_your_file.csv' with the actual path to your CSV file
df <- read.csv("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/Doctorate/Publication/Field_data/tables/Average_ARGs_per_sample.csv")

# Update "St.Donatus" to "St. Donatus" in the Location column of df
df$Location <- gsub("St.Donatus", "St. Donatus", df$Location)
# Ensure the Category column uses proper names
df$Category <- recode(df$Category, "Low" = "Low MPs", "High" = "High MPs")

# Create the boxplot with dot overlay, color by Category
boxplot_chart <- ggplot(df, aes(x=Location, y=Mean.adjCMR, fill=Category)) +
  # Boxplot for distribution
  geom_boxplot(outlier.shape = NA, alpha=0.5) +  # Hide outliers in the boxplot itself
  # Jitter to overlay individual points
  geom_jitter(aes(color=Category), width=0.2, size=2, alpha=0.7) +
  # Labels and styling
  labs(y="Average CPM per location", x="Location") +
  scale_fill_manual(values=c("Low MPs" = "#6885d0", "High MPs" = "#cb5658")) +  # Color by Category
  scale_color_manual(values=c("Low MPs" = "#6885d0", "High MPs" = "#cb5658")) +  # Color dots by Category
  theme_pubr() +  # Optionally, use theme_minimal() for a clean look
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(boxplot_chart)

# Optional: Perform a statistical test
kruskal.test(Mean.adjCMR ~ Location, data=df) # Perform the Kruskal-Wallis test

"""
Kruskal-Wallis rank sum test

data:  Mean.adjCMR by Location
Kruskal-Wallis chi-squared = 16.669, df = 7, p-value = 0.01966
"""

# Perform Dunn's post-hoc test for pairwise comparisons
dunn_test <- dunn.test(df$Mean.adjCMR, df$Location, method="BH")  # You can also use "holm" or "BH" for p-value correction
print(dunn_test) #Significance lost after multiple testing correction
"""
$P.adjusted
 [1] 0.38322213 0.46453358 0.30839969 0.49356031 0.44152617 0.48156942 0.43135422 0.26216034 0.49987311 0.45100031 0.09891246 0.30949466
[13] 0.08978270 0.25539795 0.11326417 0.07922054 0.16688497 0.15924171 0.13092611 0.23727150 0.27497981 0.48963815 0.43847523 0.45747815
[25] 0.50692420 0.43446288 0.15927548 0.09456976
"""
