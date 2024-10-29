#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 10:46:08 2024

@author: u0145079
"""

import pandas as pd
from statsmodels.multivariate.manova import MANOVA
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Load the dataset
file_path = '/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Publication/tables/Plastic-composition-Sediment.csv'  # Adjust path as necessary
data = pd.read_csv(file_path)

# Selecting only relevant "_per_gram" columns and "Category" for MANOVA
manova_data = data[['PE_PP_per_gram', 'PS_per_gram', 'PET_polyester_per_gram', 'Category']]

# Prepare the formula for MANOVA: Dependent variables ~ Independent variable
manova_formula = 'PE_PP_per_gram + PS_per_gram + PET_polyester_per_gram ~ Category'

# Run MANOVA
mv = MANOVA.from_formula(manova_formula, data=manova_data)
manova_results = mv.mv_test()

# Display the results
print(manova_results)

"""
 Multivariate linear model
============================================================
                                                            
------------------------------------------------------------
       Intercept        Value  Num DF  Den DF F Value Pr > F
------------------------------------------------------------
          Wilks' lambda 0.5905 3.0000 26.0000  6.0111 0.0030
         Pillai's trace 0.4095 3.0000 26.0000  6.0111 0.0030
 Hotelling-Lawley trace 0.6936 3.0000 26.0000  6.0111 0.0030
    Roy's greatest root 0.6936 3.0000 26.0000  6.0111 0.0030
------------------------------------------------------------
                                                            
------------------------------------------------------------
        Category        Value  Num DF  Den DF F Value Pr > F
------------------------------------------------------------
          Wilks' lambda 0.6842 3.0000 26.0000  4.0011 0.0181
         Pillai's trace 0.3158 3.0000 26.0000  4.0011 0.0181
 Hotelling-Lawley trace 0.4617 3.0000 26.0000  4.0011 0.0181
    Roy's greatest root 0.4617 3.0000 26.0000  4.0011 0.0181
============================================================
"""

#############

# Calculate mean concentrations for each microplastic type by category
mean_concentrations = data.groupby('Category')[['PE_PP_per_gram', 'PS_per_gram', 'PET_polyester_per_gram']].mean().reset_index()

# Melt the DataFrame to make it suitable for sns.barplot
melted_mean_concentrations = pd.melt(mean_concentrations, id_vars=['Category'], value_vars=['PE_PP_per_gram', 'PS_per_gram', 'PET_polyester_per_gram'],
                                     var_name='Microplastic', value_name='Mean Concentration')

# Replace the 'Microplastic' column values with the new names
replacement_names = {
    'PE_PP_per_gram': 'Polyethylene/Polypropylene',
    'PS_per_gram': 'Polystyrene',
    'PET_polyester_per_gram': 'PET/Polyester'
}
melted_mean_concentrations['Microplastic'] = melted_mean_concentrations['Microplastic'].map(replacement_names)

# Define a custom color palette
custom_palette = {'Natural pond': 'powderblue', 'Artifical pond': 'lightcoral'}

plt.figure(figsize=(8, 6))  # Adjust figure size as needed
sns.barplot(x='Microplastic', y='Mean Concentration', hue='Category', data=melted_mean_concentrations, palette=custom_palette)

# Customizing the plot's appearance
plt.title('Mean Microplastic Concentrations per Gram in Natural vs. Artificial Ponds', fontsize=14)
plt.ylabel('Mean Concentration per gram of sediment', fontsize=12)
plt.xlabel('Type of Microplastic', fontsize=12)
plt.xticks(rotation=45, fontsize=10)  # Rotate x-axis labels for better readability
plt.legend(title='Pond Category', fontsize=10)
plt.tight_layout()

# Saving the figure
plt.savefig('/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Publication/Figures/Average_MPs_composition.png', format='pdf', dpi=300)

############

# Creating a 3D scatter plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Assigning each category a color
colors = {'Natural pond': 'powderblue', 'Artifical pond': 'lightcoral'}

# Plotting each point in 3D space
for category, group in data.groupby('Category'):
    ax.scatter(group['PE_PP_per_gram'], group['PS_per_gram'], group['PET_polyester_per_gram'], 
               label=category, s=50, alpha=0.7, color=colors[category])

# Customizing axis labels with more descriptive titles
ax.set_xlabel('Polyethylene/Polypropylene')
ax.set_ylabel('Polystyrene')
ax.set_zlabel('PET/Polyester')

# Title and Legend
ax.set_title('3D Scatter Plot of Microplastics Concentrations by Pond Category')
ax.legend()

plt.show()

