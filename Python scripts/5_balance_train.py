from imblearn.over_sampling import SMOTE
from imblearn.over_sampling import ADASYN
from imblearn.under_sampling import RandomUnderSampler
from imblearn.pipeline import Pipeline
from imblearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold


import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.stats as stats
from scipy.stats import wilcoxon
from itertools import combinations
import pandas as pd
from collections import Counter


# Define strategies to balance the different training datasets

# SMOTE_1_BALANCED
def smote_1_balanced(X_train, y_train):
    ### Define resampling SMOTE_1 balanced, over and under
    #balanced_SMOTE = {2: 46, 1: 39, 3: 32, 4: 20, 5: 20, 6: 20} ### Adjust to use oversample
    over = SMOTE(random_state=42, sampling_strategy= {4: 20, 5: 20, 6: 20}, k_neighbors=3)
    #balance_under = {2: 25, 1: 25, 3: 25, 4: 20, 5: 20, 6: 20} ### Adjust to use undersample
    under = RandomUnderSampler(random_state=42, sampling_strategy={2: 25, 1: 25, 3: 25})
    steps = [('o', over), ('u', under)]
    pipeline = Pipeline(steps=steps)
    X_train_smote, y_train_smote = pipeline.fit_resample(X_train, y_train)
    return X_train_smote, y_train_smote

# SMOTE_2_BALANCED
def smote_2_balanced(X_train, y_train):
    ### Define resampling  SMOTE_2
    #balanced_SMOTE = {1: 15, 2: 15, 3: 15, 4: 10, 5: 10, 6: 10} ### Adjust to use oversample
    over = SMOTE(random_state=42, sampling_strategy= {4: 15, 6: 15}, k_neighbors=3)
    #balance_under = = {1: 15, 2: 15, 3: 15, 4: 10, 5: 10, 6: 10} ### Adjust to use undersample
    under = RandomUnderSampler(random_state=42, sampling_strategy={1: 15, 2: 15, 3: 15, 5: 15})
    steps = [('o', over), ('u', under)]
    pipeline = Pipeline(steps=steps)
    X_train_smote, y_train_smote = pipeline.fit_resample(X_train, y_train)
    return X_train_smote, y_train_smote
    return X_smote2, y_smote2

# ADASYN_1_BALANCED
def adasyn_1_balanced(X_train, y_train):
    ### Define resampling ADASYN_1 balanced, over and under
    #balanced_Adasyn = {2: 46, 1: 39, 3: 32, 4: 20, 5: 20, 6: 20} ### Adjust to use oversample
    over = ADASYN(random_state=42, sampling_strategy= {4: 20, 5: 20, 6: 20}, n_neighbors=4)
    #balance_under = {2: 25, 1: 25, 3: 25, 4: 20, 5: 20, 6: 20} ### Adjust to use undersample
    under = RandomUnderSampler(random_state=42, sampling_strategy={2: 25, 1: 25, 3: 25})
    steps = [('o', over), ('u', under)]
    pipeline = Pipeline(steps=steps)
    X_train_smote, y_train_smote = pipeline.fit_resample(X_train, y_train)
    return X_train_smote, y_train_smote

# ADASYN_2_BALANCED
def adasyn_2_balanced(X_train, y_train):
    ### Define resampling default ADASYN_2
    #balanced_Adasyn = {2: 15, 1: 15, 3: 15, 4: 10, 5: 10, 6: 10} ### Adjust to use oversample
    over = ADASYN(random_state=42, sampling_strategy= {4: 15, 6: 15}, n_neighbors=4)
    #balance_under = {2: 15, 1: 15, 3: 15, 4: 10, 5: 10, 6: 10} ### Adjust to use undersample
    under = RandomUnderSampler(random_state=42, sampling_strategy={1: 15, 2: 15, 3: 15, 5: 15})
    ### By default ADASYN equalizes all the classes to the same number as the majority class
    steps = [('o', over), ('u', under)]
    pipeline = Pipeline(steps=steps)
    ### random_state is the seed used by the random number generator, k_neighbors is the number of nearest neighbors to used to construct synthetic samples
    X_ada2, y_ada2 = pipeline.fit_resample(X_train, y_train)
    return X_ada2, y_ada2

## Function to plot_histogram 
def plot_histogram(ax, Patient_subgroups, title):
    """
    Create a histogram of the gene expression data for a list of genes
    
    Parameters:
    ax (matplotlib.axes._subplots.AxesSubplot): The axis where the histogram will be plotted.
    data (Patient subgroups): Dictionary with the patient subgroups, keys patients and values # of patients in that subgroup
    title (str): The title of the subplot
    
    Returns:
    None
    """
    counter = Counter(Patient_subgroups)
    for k, v in counter.items():
        per = v / len(Patient_subgroups) * 100
        print('Class=%d, n=%d (%.3f%%)' % (k, v, per))
    
    # Plot the distribution using bar chart on the current axis
    ax.bar(counter.keys(), counter.values())
    
    # Customize the plot with title and axis labels
    ax.set_title(title)
    ax.set_xlabel('Patient subgroups')
    ax.set_ylabel('Number of patients')

    
    
def get_metrics(results_dict):
    # Create an empty dataframe
    df = pd.DataFrame()
    # Iterate through the results_dict
    for model_name, scores in results_dict.items():
        for score in scores:
            #Create a temporary datafram
            temp_df = pd.DataFrame.from_dict({
                'model': model_name,
                'test_Matthews_corrcoef': score['test_Matthews_corrcoef'],
                'test_balanced_accuracy': score['test_balanced_accuracy'],
                'test_f1_score': score['test_f1_score'],
                'test_precision': score['test_precision'],
                'test_recall': score['test_recall'],
                'fit_time': score['fit_time'],
                'score_time': score['score_time']
            })
            # Append it to the main dataframe
            df = df.append(temp_df, ignore_index=True)
    return df

def set_names(df):
    dataset_names = {0: 'Imbalanced', 1: 'SMOTE_1', 2: 'SMOTE_2', 3: 'ADASYN_1', 4: 'ADASYN_2'}
    model_names = {'RandomForestClassifier': 'RF','DecisionTreeClassifier': 'DT','XGBClassifier': 'XGBoost','SVC': 'SVC','KNeighborsClassifier': 'KNN'}
    df['dataset'].replace(dataset_names, inplace=True)
    df['model'].replace(model_names, inplace=True)
    return df


def plot_boxplot_metrics(results_all):
    '''Plot Matthews corrcoef and balanced accuracy scores'''
    ### Set figure size
    plt.figure(figsize=(20, 10))

    # Matthews corrcoef
    plt.subplot(1, 2, 1)
    sns.boxplot(x='model', y='test_Matthews_corrcoef', hue='dataset', data=results_all)
    plt.xlabel('Model')
    plt.ylabel('MCC')
    plt.title('Matthews Corrcoef scores')
    plt.xticks(rotation=90)
    
    # Save the plot as a high-resolution image (e.g., PNG or PDF)
    plt.savefig('results_metrics.png', dpi=300, bbox_inches='tight')




def calculate_mean_std(df, group_cols, value_cols, title=None):
    grouped_df = df.groupby(group_cols).agg({col: ['mean', 'std'] for col in value_cols})
    grouped_df.reset_index(inplace=True)
    # Save the grouped DataFrame to a CSV file if a title is provided
    if title:
        grouped_df.to_csv(title, index=False)  # Set index=False to exclude the index column
    return grouped_df


def bonferroni_correction(p_values):
    num_tests = len(p_values)
    corrected_p_values = [p * num_tests for p in p_values]
    return corrected_p_values



def calculate_validation_stats(results_all):
    
    # Extract the necessary columns
    df2 = results_all[['dataset','test_Matthews_corrcoef']]
    
    # Perform Wilcoxon pairwise comparisons
    groups = df2['dataset'].unique()
    pairs = list(combinations(groups, 2))
    
    wilcoxon_data = []
    for group1, group2 in pairs:
        data1 = df2[df2['dataset'] == group1]['test_Matthews_corrcoef']
        data2 = df2[df2['dataset'] == group2]['test_Matthews_corrcoef']
        statistic, p_value = wilcoxon(data1, data2, zero_method='pratt')
        wilcoxon_data.append([group1, group2, statistic, p_value])
    
    # Create DataFrame for the Wilcoxon results
    wilcoxon_df = pd.DataFrame(wilcoxon_data, columns=['Dataset 1', 'Dataset 2', 'Wilcoxon Statistic', 'Wilcoxon p-value'])
    
    # Apply Bonferroni correction to p-values
    wilcoxon_df['Corrected p-value'] = bonferroni_correction(wilcoxon_df['Wilcoxon p-value'])
    
    return wilcoxon_df



def plot_boxplot_metrics(results_all):
    # Set the font size and font family
    plt.rcParams.update({'font.size': 12, 'font.family': 'serif'})

    # Create a figure with the specified size and aspect ratio
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Matthews corrcoef
    sns.boxplot(x='model', y='test_Matthews_corrcoef', hue='dataset', data=results_all, palette='Set3', ax=axes[0])
    axes[0].set_xlabel('Model')
    axes[0].set_ylabel('MCC')
    #axes[0].set_title('Matthews correlation coefficient')
    axes[0].tick_params(axis='x', rotation=45)
    axes[0].grid(True, linestyle='--', alpha=0.7)
    axes[0].set_ylim(0,1)
    axes[0].get_legend().remove()

    # Balanced accuracy
    sns.boxplot(x='model', y='test_balanced_accuracy', hue='dataset', data=results_all, palette='Set3', ax=axes[1])
    #axes[1].set_title('Balanced Accuracy')
    axes[1].set_xlabel('Model')
    axes[1].set_ylabel('BA')
    axes[1].tick_params(axis='x', rotation=45)
    axes[1].grid(True, linestyle='--', alpha=0.7)
    axes[1].set_ylim(0,1)

    # Adjust spacing between subplots
    plt.tight_layout(pad=2.0)
    
    # Add labels "A" and "B" to each subplot
    axes[0].text(-0.15, 1.1, "F", transform=axes[0].transAxes, fontsize=14, fontweight='bold')
    axes[1].text(-0.15, 1.1, "G", transform=axes[1].transAxes, fontsize=14, fontweight='bold')


    # Move the legend outside and to the upper right corner
    plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')

    
         
    # Save the plot as a high-resolution image (e.g., PNG or PDF)
    plt.savefig('results_plot.png', dpi=800, bbox_inches='tight')
    
    # Show the plot (optional)
    plt.show()

# Example usage:
# plot_boxplot_metrics(results_all)



