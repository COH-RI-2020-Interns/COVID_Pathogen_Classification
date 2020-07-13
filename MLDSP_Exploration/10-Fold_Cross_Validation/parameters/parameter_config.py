import json
from os import getcwd, listdir

parameter_dict = {}

# Logistic Regression Parameters

lr = {'penalty': ['l1', 'l2', 'elasticnet', None],
        'C': [0.001, 0.005, 0.2, 0.6, 1],
        'solver': 'liblinear',
        'class_weight': [None, 'balanced'],
        }

parameter_dict['log'] = lr

# Random Forest Parameters

rf = {'n_estimators': [100, 150, 200, 350, 500],
      'criterion': ['gini', 'entropy'],
      'max_depth': [None, 3, 5]
        }

parameter_dict['rf'] = rf

# AdaBoost Parameters

ada = {'n_estimators': [50, 75, 100]}

parameter_dict['ada'] = ada

# K Neighbors Parameters

knn = {'n_neighbors': [5, 7, 9, 11],
        'weights': ['uniform', 'distance'],
        'p': [1, 2]
        }

parameter_dict['knn'] = knn

# Support Vector Machine Parameters

svm = {'C': [0.001, 0.005, 0.2, 0.6, 1],
       'kernel': ['linear', 'rbf', 'polynomial'],
       'class_weight': [None, 'balanced']
       }

parameter_dict['svm'] = svm

# Decision Tree Classifier Parameter

decision_tree = {'criterion': ['gini', 'entropy'], #default = 'gini'
        'splitter': ['best', 'random'], #default = 'best'
        'class_weight': [None]}

parameter_dict['decision_tree'] = decision_tree

data_path = getcwd() + "/data"

with open(f"{data_path}/{sorted(listdir(data_path))[0][2]}", "w") as my_file:
    json.dump(parameter_dict, my_file)
