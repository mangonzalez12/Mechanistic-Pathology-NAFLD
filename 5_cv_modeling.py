import pandas as pd
from sklearn.model_selection import KFold, RandomizedSearchCV, cross_validate, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from xgboost import XGBClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from modeling_preprocessing import filter_df, get_metrics
from sklearn.metrics import make_scorer, matthews_corrcoef, balanced_accuracy_score, precision_score, recall_score, f1_score
from sklearn.model_selection import StratifiedKFold

import imblearn
from imblearn.over_sampling import SMOTE
from imblearn.over_sampling import ADASYN
from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import RandomUnderSampler
from imblearn.pipeline import make_pipeline
from imblearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split
import warnings
import scipy.stats as stats
from scipy.stats import wilcoxon
from itertools import combinations

# Suppress the specific warning
warnings.filterwarnings("ignore", message="Exact p-value calculation does not work", category=UserWarning)

from balance_train import smote_1_balanced, smote_2_balanced, adasyn_1_balanced, adasyn_2_balanced, set_names
from balance_train import plot_boxplot_metrics, calculate_mean_std, bonferroni_correction, calculate_validation_stats

#Read in data
### Read input data and transpose to have genes as columns, note, last column is the label for each datapoint ("Cluster WGCNA")
df = pd.read_csv("Govaere_2rlog_visualization_6_subgroups WGCNA_800genes.csv", index_col=0)
#label = df.iloc[:,-1]

##merge label with df
#df = df.merge(label, left_index=True, right_index=True)

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


# Define the datasets
datasets = [
    {'X_train': X_train_1, 'y_train': y_train_1},
    {'X_train': X_train_2, 'y_train': y_train_2},
    {'X_train': X_train_3, 'y_train': y_train_3},
    {'X_train': X_train_4, 'y_train': y_train_4},
    {'X_train': X_train_5, 'y_train': y_train_5}
]

# Define the models and their hyperparameters
models = {
    'RandomForestClassifier': RandomForestClassifier(),
    'DecisionTreeClassifier': DecisionTreeClassifier(),
    'XGBClassifier': XGBClassifier(),
    'KNeighborsClassifier': KNeighborsClassifier()
}

hyperparameters = {
    'RandomForestClassifier': {
        'max_depth': [4, 6, 8, 10],
        'n_estimators': [200, 400, 600, 800, 1000],
        'bootstrap': [True, False],
        'max_features': ['log2', 'sqrt']
    },
    'DecisionTreeClassifier': {
        'max_depth': [2,3,4,5],
        'max_features': ['sqrt', 'log2'],
        'min_samples_split': [2,3,4,5],
        'min_samples_leaf': [1,2,3,4,5],
        'max_leaf_nodes': [2,3,4,5,6,7,8,9,10],
        'criterion': ['gini', 'entropy']
    },
    'XGBClassifier': {
        'max_depth': [4, 6, 8, 10],
        'n_estimators': [200, 400, 600, 800, 1000],
        "colsample_bytree": [0.3, 0.7, 1.0],
        "subsample": [0.2, 0.4, 0.6, 0.8, 1],
        "min_child_weigth": [0.6, 0.8],
        "min_split_loss": [0.2, 0.4], 
        "num_class": [6]#Number of classes in the dataset
    },
    'KNeighborsClassifier': {
        'n_neighbors': [3,4,5,6],
        'weights': ['uniform', 'distance'],
        'algorithm': ['auto', 'ball_tree', 'kd_tree']
    }
}

n_jobs = 30
###Declare the scoring metrics
scoring = {'Matthews_corrcoef': make_scorer(matthews_corrcoef),
           'balanced_accuracy': make_scorer(balanced_accuracy_score),
           'f1_score': make_scorer(f1_score, average='micro'),
           'precision': make_scorer(precision_score, average='micro'),
           'recall': make_scorer(recall_score, average='micro')       
}

###Create a dataframe to store all the results
results_all = pd.DataFrame()


# Loop through each dataset
for idx, data in enumerate(datasets):
    # Define the input data
    X_train = data['X_train']
    y_train = data['y_train']


    ###Create a dictionary to store the results 
    results_dict = {}  ##

    score_names = ['Matthews_corrcoef', 'balanced_accuracy', 'f1_score', 'precision', 'recall']


    ###Fitting and evaluating each model
    for model_name, model in models.items():
        ###Declare the inner and outer cross validation
        inner_cv =  StratifiedKFold(n_splits=2, shuffle=True, random_state=1)
        outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)
    
        random_search = RandomizedSearchCV(model, hyperparameters[model_name], cv=inner_cv, n_jobs=n_jobs)
        for i in range(10):
            scores = cross_validate(random_search, X_train, y_train, cv=outer_cv, scoring=scoring, n_jobs=n_jobs)
            ###Store the results in the dictionary
            #model_name = ['RandomForestClassifier','DecisionTreeClassifier','XGBClassifier','KNeighborsClassifier', 'SVC']
            #score_names = ['CV Matthews_corrcoef', 'CV balanced_accuracy',  'CV precision', 'CV recall']
            print(f"Round {i+1} scores: {scores}")
            # Store the results in the results_dict
            if model_name not in results_dict:
                results_dict[model_name] = []
            results_dict[model_name].append(scores)
            
    results_df = get_metrics(results_dict)
    ##Add a column to the dataframe to indicate the dataset
    results_df['dataset'] = idx
    

    
    ##Append the results to the combined dataframe
    results_all = results_all.append(results_df, ignore_index=True)

results_all = set_names(results_all)

###Print the combined results
results_all.to_csv("results_multiclass.csv", index=False)



# Plot results
Plot_6 = plot_boxplot_metrics(results_all)

# Calculate mean and standard deviation
# Group values
value_cols = ['test_Matthews_corrcoef', 'test_balanced_accuracy', 'test_f1_score', 'test_precision', 'test_recall']
group_cols=['model', 'dataset']
# Write dataframe
stats_df = calculate_mean_std(df=results_all, group_cols=group_cols, value_cols=value_cols)
stats_df.to_csv("Mean and standard deviation Results.csv")


