import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List, Dict
import pandas as pd
from sklearn.model_selection import StratifiedKFold, RandomizedSearchCV, cross_validate

from sklearn.manifold import TSNE
from matplotlib.patches import Patch
from collections import Counter

from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from xgboost import XGBClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.metrics import make_scorer, matthews_corrcoef, balanced_accuracy_score, f1_score, precision_score, recall_score
from sklearn.model_selection import cross_validate
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV

from imblearn.over_sampling import SMOTE
from imblearn.over_sampling import ADASYN
from imblearn.under_sampling import RandomUnderSampler
from imblearn.pipeline import Pipeline
from imblearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold


## Preprocessing before modeling for feature selection

##Function to filter datasets on particular gene list
def filter_df(df, gene_list):
    """
    Filter a dataframe of gene expression data with a list of genes
    
    Parameters:
    df (DataFrame): The gene expression data
    gene_list (list): The list of genes to be used for filtering
    
    Returns:
    filtered_df (DataFrame): The filtered dataframe
    """
    
    # Filter the dataframe columns based on the gene list
    filtered_df = df.filter(items=gene_list)
    
    # Return the filtered dataframe
    return filtered_df


def plot_scatterplot(x, y, xlabel, ylabel, title):
    sns.set_style("darkgrid") 
    sns.scatterplot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()


##Function to remove correlated features
def remove_correlated_features(dataframe, threshold):
    col_corr = set() # Set of all the names of correlated columns
    corr_matrix = dataframe.corr()
    for i in range(len(corr_matrix.columns)):
        for j in range(i):
            if abs(corr_matrix.iloc[i, j]) > threshold:
                # getting the name of column
                colname = corr_matrix.columns[i]
                col_corr.add(colname)
    filtered_dataframe = dataframe.drop(columns = col_corr)
    return filtered_dataframe


## Function to create a heatmap of certain genes
def heatmap(dataframe, metadata, gene_list_df, title):
    """
    Create a heatmap of the gene expression data for a list of genes
    
    Parameters:
    dataframe (DataFrame): The gene expression data, genes in columns
    gene_list (list): The list of genes to be used for filtering
    
    Returns:
    None
    """
    # Filter the dataframe columns based on the gene list
    meta_genes = gene_list_df
    gene_list = meta_genes['Ens_id'].tolist()
    
    ###Rownames for the heatmap are the gene names
    #meta_genes = gene_list_df['Hub_gene']
    
    filtered_df = dataframe.filter(items=gene_list)
    trns_df = filtered_df.T
    # Add the metadata to the filtered dataframe to set the colors of Fibrosis scores   
    Fibrosis_labels = metadata['Fibrosis_label']
    lut = dict(Fibrosis_0= "#FADDC3",Fibrosis_1= "#F9C29C", Fibrosis_2= "#F6A173", Fibrosis_3= "#F17B51", Fibrosis_4= "#EA4C3B")
    Fibrosis_colors = pd.Series(Fibrosis_labels).map(lut)
    
    # Add the metadata to the filtered dataframe to set the colors of nas_score
    Nas_scores = metadata['nas_score']
    lut2 = dict(nas_score_0 = "#F4FAFE", nas_score_1 = "#DEEEF7", nas_score_2 = "#C1DBEC", nas_score_3 = "#A1C4E0", nas_score_4 = "#7FABD3", nas_score_5 = "#5C90C6" , nas_score_6 = "#3573B9",  nas_score_7 = "#305596",  nas_score_8 = "#273871")  
    Nas_colors = pd.Series(Nas_scores).map(lut2)
    
    column_cols = pd.concat([Fibrosis_colors, Nas_colors], axis=1)
    
    # Create a heatmap using clustermap, metric 'Euclidean' on both cols and rows of the filtered dataframe
    fig = sns.clustermap(trns_df, col_colors=column_cols, cmap='coolwarm', standard_scale=None, metric='euclidean',
                         col_cluster=True, row_cluster=True)## Cluster on cols and rows
    
    ### Add the legends to the heatmap
    leg_fib = [Patch(facecolor=lut[Fibrosis_colors]) for Fibrosis_colors in lut]
    legend = fig.ax_heatmap.legend(leg_fib, lut.keys(), title="Fibrosis_label", bbox_to_anchor=(1.5, 1.3), loc=1)
    leg_nas = [Patch(facecolor=lut2[Nas_colors]) for Nas_colors in lut2]
    fig.ax_col_colors.legend(leg_nas, lut2.keys(), title="nas_score", bbox_to_anchor=(1.6, 0.5), loc=3)

    ###Save the heatmap
    fig.savefig('heatmap.png', dpi=500)
    plt.title(title)
    fig.figsize=(40, 20)
    plt.xlabel('Genes')
    plt.show()
    
    
## Function to create a heatmap of certain genes
def plot_histogram(Patient_subgroups, title):
    """
    Create a histogram of the gene expression data for a list of genes
    
    Parameters:
    data (Patient subgroups): Dictionary with the patient subgroups, keys patients and values # of patients in that subgroup
    
    Returns:
    None
    """
    counter = Counter(Patient_subgroups)
    for k,v in counter.items():
        per = v / len(Patient_subgroups) * 100
        print('Class=%d, n=%d (%.3f%%)' % (k, v, per))
    # plot the distribution
    plt.bar(counter.keys(), counter.values())
    plt.title(title)
    plt.xlabel('Patient subgroups')
    plt.ylabel('Number of patients')
    plt.show()
    
def evaluate_model(y_test, y_pred):
    acc_score =  accuracy_score(y_test, y_pred)
    cf_mat = confusion_matrix(y_test, y_pred)
    classification_report =  classification_report(y_test, y_pred)
    return acc_score, cf_mat, classification_report
    
    
def plot_tsne(dataframe, metadata, title):
    """
    Create a t-SNE plot of the gene expression data
    
    Parameters:
    dataframe (DataFrame): The gene expression data, genes (dimensions) in columns to reduce to 2 dimensions
    
    Returns:
    None
    """
    # Create a t-SNE plot of the dataframe
    tsne = TSNE(n_components=2, perplexity=30, learning_rate=200, n_iter=1000, random_state = 0)##Default parameters
    tsne_df = tsne.fit_transform(dataframe)
    tsne_df = pd.DataFrame(tsne_df, columns=['x', 'y'])
    tsne_df['label'] = metadata['Fibrosis_label']
    
    # Create a scatterplot of the t-SNE plot
    sns.scatterplot(x='x', y='y', hue='label', data=tsne_df, palette='coolwarm', legend='full')
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.title(title)
    plt.show()


def smote_1_balanced(X_train, y_train):
    ### Define resampling SMOTE_1 balanced, over and under
    #balanced_SMOTE = {2: 46, 1: 39, 3: 32, 4: 20, 5: 20, 6: 20} ### Adjust to use oversample
    over = SMOTE(random_state=42, sampling_strategy= {4: 20, 5: 20, 6: 20})
    #balance_under = {2: 25, 1: 25, 3: 25, 4: 20, 5: 20, 6: 20} ### Adjust to use undersample
    under = RandomUnderSampler(random_state=42, sampling_strategy={2: 25, 1: 25, 3: 25})
    steps = [('o', over), ('u', under)]
    pipeline = Pipeline(steps=steps)
    X_train_smote, y_train_smote = pipeline.fit_resample(X_train, y_train)
    return X_train_smote, y_train_smote

def smote_2_default(X_train, y_train):
    ### Define resampling default SMOTE_2
    ### By default SMOTE equalizes all the classes to the same number as the majority class, if strategy not specified
    ### Adjust to majority class
    sm = SMOTE(random_state=42) 
    ### random_state is the seed used by the random number generator, k_neighbors is the number of nearest neighbors to used to construct synthetic samples
    X_smote2, y_smote2 = sm.fit_resample(X_train, y_train)
    return X_smote2, y_smote2

def adasyn_1_balanced(X_train, y_train):
    ### Define resampling ADASYN_1 balanced, over and under
    #balanced_Adasyn = {2: 46, 1: 39, 3: 32, 4: 20, 5: 20, 6: 20} ### Adjust to use oversample
    over = ADASYN(random_state=42, sampling_strategy= {4: 20, 5: 20, 6: 20})
    #balance_under = {2: 25, 1: 25, 3: 25, 4: 20, 5: 20, 6: 20} ### Adjust to use undersample
    under = RandomUnderSampler(random_state=42, sampling_strategy={2: 25, 1: 25, 3: 25})
    steps = [('o', over), ('u', under)]
    pipeline = Pipeline(steps=steps)
    X_train_smote, y_train_smote = pipeline.fit_resample(X_train, y_train)
    return X_train_smote, y_train_smote

def adasyn_2_default(X_train, y_train):
    ### Define resampling default ADASYN_2
    ### By default ADASYN equalizes all the classes to the same number as the majority class
    ### Adjust to majority class
    ada_sm = ADASYN(random_state=42, sampling_strategy='all')
    ### random_state is the seed used by the random number generator, k_neighbors is the number of nearest neighbors to used to construct synthetic samples
    X_ada2, y_ada2 = ada_sm.fit_resample(X_train, y_train)
    return X_ada2, y_ada2


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

    # Balanced accuracy
    plt.subplot(1, 2, 2)
    sns.boxplot(x='model', y='test_balanced_accuracy', hue='dataset', data=results_all)
    plt.title('Balanced Accuracy scores')
    plt.xlabel('Model')
    plt.ylabel('Balanced Accuracy')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()
    
    def plot_boxplot_metrics2(results_all):
        '''Plot F1 and Precision scores'''
    ### Set figure size
    plt.figure(figsize=(20, 10))

    # F1 score
    plt.subplot(1, 2, 1)
    sns.boxplot(x='model', y='test_f1_score', hue='dataset', data=results_all)
    plt.xlabel('Model')
    plt.ylabel('F1 score')
    plt.title(' F1 scores')
    plt.xticks(rotation=90)

    # Precision
    plt.subplot(1, 2, 2)
    sns.boxplot(x='model', y='test_precision', hue='dataset', data=results_all)
    plt.title('Precision scores')
    plt.xlabel('Model')
    plt.ylabel('Precision')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()



def run_model(datasets: List[Dict[str, pd.DataFrame]], models: Dict[str, object], hyperparameters: Dict[str, Dict[str, object]], n_jobs: int = 30) -> pd.DataFrame:
    scoring = {'Matthews_corrcoef': make_scorer(matthews_corrcoef),
           'balanced_accuracy': make_scorer(balanced_accuracy_score),
           'f1_score': make_scorer(f1_score, average='micro'),
           'precision': make_scorer(precision_score, average='micro'),
           'recall': make_scorer(recall_score, average='micro')       
    }
    results_all = pd.DataFrame()

    # Loop through each dataset
    for idx, data in enumerate(datasets):
        # Define the input data
        X_train = data['X_train']
        y_train = data['y_train']

        # Create a dictionary to store the results 
        results_dict = {}

        score_names = ['Matthews_corrcoef', 'balanced_accuracy', 'f1_score', 'precision', 'recall']

        # Fitting and evaluating each model
        for model_name, model in models.items():
            # Declare the inner and outer cross validation
            inner_cv = StratifiedKFold(n_splits=2, shuffle=True, random_state=1)
            outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)

            # Perform hyperparameter tuning using randomized search
            random_search = RandomizedSearchCV(model, hyperparameters[model_name], cv=inner_cv, n_jobs=30)

            # Fit and evaluate the model with cross-validation
            for i in range(1):
                scores = cross_validate(random_search, X_train, y_train, cv=outer_cv, scoring=scoring)

                # Store the results in the dictionary
                if model_name not in results_dict:
                    results_dict[model_name] = []
                results_dict[model_name].append(scores)

        # Convert the results dictionary to a dataframe
        results_df = pd.DataFrame()

        for model_name, results_list in results_dict.items():
            for i, result in enumerate(results_list):
                for score_name, score in result.items():
                    results_df.loc[f'{model_name}_{i}', score_name] = score.mean()

        # Add a column to the dataframe to indicate the dataset
        results_df['dataset'] = idx

        # Append the results to the combined dataframe
        results_all = results_all.append(results_df, ignore_index=True)
        ###Print the combined results
        results_all

    return results_all.to_csv("results_all_genes_combat_multiclass.csv", index=False)




    
def plot_scatterplot_DEGs(x, y, DEGs_df, title, threshold_positive, threshold_negative):
        '''Plot scatterplot of DEGs with annotated genes based on the log2foldchange especified in the thresholds'''
        ### Set figure size
        plt.figure(figsize=(20, 10))
        ### Plot genes with log2foldchange higher than 2 and lower than 2
        fig, axs = plt.subplots(ncols=2, figsize=(20, 10))
        sns.scatterplot(x=x, y=y, data=DEGs_df, palette="colorblind", ax=axs[0])
        for i, txt in enumerate(Gene.name):
            if x[i] > threshold_positive or x[i] < threshold_negative and y[i] > threshold_positive or y[i] < threshold_negative:
                axs[0].annotate(txt, (x[i], y[i]))


