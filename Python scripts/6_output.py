import scipy.stats as stats
from scipy.stats import wilcoxon
from itertools import combinations
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import warnings
import os

# Import functions from balance_train
from balance_train import plot_boxplot_metrics, calculate_mean_std, bonferroni_correction, calculate_validation_stats
from balance_train import smote_1_balanced, smote_2_balanced, adasyn_1_balanced, adasyn_2_balanced, plot_histogram

# Suppress the specific warning
warnings.filterwarnings("ignore", message="Exact p-value calculation does not work", category=UserWarning)




# Read results file
results_all = pd.read_csv("results_all_genes_800genes_multiclass.csv")

# Plot results
Plot_6 = plot_boxplot_metrics(results_all)

# Calculate mean and standard deviation
# Group values
value_cols = ['test_Matthews_corrcoef', 'test_balanced_accuracy', 'test_f1_score', 'test_precision', 'test_recall']
group_cols=['model', 'dataset']
# Write dataframe
stats_df = calculate_mean_std(df=results_all, group_cols=group_cols, value_cols=value_cols)
stats_df.to_csv("Mean and standard deviation Results.csv")


# Choose model to calculate validation stats
df2 =results_all.loc[results_all['model']=='RF']
df3 = calculate_validation_stats(results_all=df2)
df3.to_csv("Validationn stats results.csv")


# Write files with the data balancing techniques
df = pd.read_csv("Govaere_2rlog_visualization_6_subgroups WGCNA_800genes.csv", index_col=0)

### Set up X and y
y = df["Cluster_WGCNA"]
X = df.drop("Cluster_WGCNA", axis=1)

### Data split into train and test set as follows using train_test_split function
##Used 30% of the data for testing, to allow for more data for SMOTE and test
X_train, X_test, y_train, y_test = train_test_split(X,y, random_state=0, test_size=0.30) 


####Split the data into 5 datasets
X_train_1, y_train_1 = X_train, y_train  ##imbalanced dataset
X_train_2, y_train_2 = smote_1_balanced(X_train, y_train) ##SMOTE 1
X_train_3, y_train_3 = smote_2_balanced(X_train, y_train) ##SMOTE 2
X_train_4, y_train_4 = adasyn_1_balanced(X_train, y_train) ##ADASYN 1
X_train_5, y_train_5 = adasyn_2_balanced(X_train, y_train) ##ADASYN 2




## Plot all graphs combined
# Create a 2x3 grid of subplots
fig, axs = plt.subplots(2, 3, figsize=(12, 8))

# Call your custom function for plotting histograms on each subplot
plot_histogram(axs[0,0], y_train_1, "Imbalanced")
plot_histogram(axs[0,1], y_train_2, "SMOTE 1")
plot_histogram(axs[0,2], y_train_3, "SMOTE 2")
plot_histogram(axs[1,0], y_train_4, "ADASYN 1")
plot_histogram(axs[1,1], y_train_5, "ADASYN 2")

# Hide non-existent plot
fig.delaxes(axs[1,2])

# Add labels to subplots using ax.annotate
for i, label in enumerate(['A', 'B', 'C', 'D', 'E']):
    axs.flatten()[i].annotate(label, xy=(-0.05, 1.05), xycoords='axes fraction', fontsize=12, weight='bold', va='center', ha='center')


# Adjust layout
plt.tight_layout()

# Save figure
plt.savefig('figure_subplots.png', dpi=800, bbox_inches='tight')

# Show plot
plt.show()
