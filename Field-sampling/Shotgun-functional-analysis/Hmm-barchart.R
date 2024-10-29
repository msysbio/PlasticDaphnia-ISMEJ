library(ggplot2)
library(dplyr)
library(tidyr)
library(broom) # for tidy()

# Load the data
df <- read.csv("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Publication/Field_data/tables/Enzyme-abundances-samples.csv",row.names=1)
meta = df %>%
  mutate("Alpha.beta.hydrolase.family"=NULL,"Esterase.PHB.depolymerase"=NULL,"Lipase"=NULL)

meta = meta %>%
  mutate("Location"= NULL) #REMOVE LOCATION

df = df %>%
  mutate("Location" =NULL, "Category" =NULL)

# Transform data from wide to long format
long_data <- gather(df, Key, Value, -Category)

# Calculate average values for each enzyme family by pond type
average_values <- long_data %>%
  group_by(Category, Key) %>%
  summarize(Average = mean(Value))

average_values$Key <- factor(average_values$Key, levels = unique(average_values$Key), labels = c("Alpha-beta hydrolase", "Carboxylesterase", "Lipase"))

# Now plot with the updated 'Key' column
bar_chart <- ggplot(average_values, aes(x=Key, y=Average, fill=Category)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c("blue", "red")) + # Specify your colors here
  labs(y="Average Value", fill="Pond Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), # Adjust text angle and justification for readability
        axis.text.y = element_text(),
        legend.title = element_text(face = "bold"))+
  theme_pubr()
bar_chart
# Perform the Mann-Whitney U test for each enzyme family
#Alpha-beta hydrolase
shapiro.test(df$Alpha.beta.hydrolase.family[df$Category=="High"])
shapiro.test(df$Alpha.beta.hydrolase.family[df$Category=="Low"])

wilcox.test(Alpha.beta.hydrolase.family ~ Category,data=df) #>0.05

#Carboxylesterase
shapiro.test(df$Esterase.PHB.depolymerase[df$Category=="High"])
shapiro.test(df$Esterase.PHB.depolymerase[df$Category=="Low"])

wilcox.test(Esterase.PHB.depolymerase ~ Category,data=df) #>0.05

#Lipase
shapiro.test(df$Lipase[df$Category=="High"])
shapiro.test(df$Lipase[df$Category=="Low"])

wilcox.test(Lipase ~ Category,data=df) #>0.05



