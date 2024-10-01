# Import necessary libraries
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from keras.models import Sequential
from keras.layers import Dense
from sklearn.model_selection import StratifiedKFold
from keras.wrappers.scikit_learn import KerasClassifier
import xgboost as xgb
import copy

# Function to run the simulation for a given iteration
def run_ps(i, conf_list, cov_list, data_template, output_file, random_seed=43):
    """
    Function to run propensity scores simulation on covariates, confounders, and treatment models
    
    Parameters:
    i (int): Simulation index
    conf_list (str): Path to confounder list
    cov_list (str): Path to covariate list
    data_template (str): Path template to load the data for each iteration
    output_file (str): Path template to save the output CSV
    random_seed (int): Random seed for reproducibility
    """
    
    # Load dataset for the current simulation
    simdata = pd.read_csv(data_template.format(i), index_col=0).dropna()

    # Create a copy of the simulation data
    data = pd.DataFrame(copy.deepcopy(simdata))

    # Filter confounders and covariates based on the provided lists
    conf_cols = data.filter(map(str, conf_list['x'].to_list()))
    cov_cols = data.filter(map(str, cov_list['x'].to_list()))

    # Define treatment variable
    t = data["EXPOSURE{}".format(i)]

    # Split the data (50% train, 50% test), stratified by treatment
    X_train, X_test, y_train, y_test = train_test_split(cov_cols, t, test_size=0.5, stratify=t, random_state=random_seed)

    # Fit reference logistic regression model (with confounders only)
    logmod_con = LogisticRegression(penalty='none', solver='lbfgs', max_iter=1000)
    logmod_con.fit(conf_cols, t)

    # Cross-validation setup (10 folds)
    skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=1001)

    # ------------------ LASSO Logistic Regression ------------------ #
    logistic = LogisticRegression(penalty='l1', max_iter=1000)
    alphas = np.arange(0.01, 0.31, 0.02)
    hyperparameters = {'C': alphas, 'solver': ['liblinear', 'saga']}

    # Randomized search for hyperparameter tuning with LASSO
    clf_ps = RandomizedSearchCV(estimator=logistic, param_distributions=hyperparameters, 
                                cv=skf.split(X_train, y_train), verbose=0,
                                scoring='neg_brier_score', n_iter=100)
    logmod_ps = clf_ps.fit(X_train, y_train)

    # Fit the tuned model on full covariate data
    logmod_ps.fit(cov_cols, t)

    # ------------------ Multi-Layer Perceptron (MLP) ------------------ #
    def build_classifier(optimizer, kernel, units, hidden_layers, activation):
        classifier = Sequential()
        # Input and first hidden layer
        classifier.add(Dense(units=units, activation='sigmoid', input_dim=cov_cols.shape[1], kernel_initializer=kernel))
        # Adding more hidden layers
        for _ in range(hidden_layers):
            classifier.add(Dense(units=units, activation=activation, kernel_initializer=kernel))
        # Output layer (binary classification)
        classifier.add(Dense(1, activation='sigmoid', kernel_initializer=kernel))
        classifier.compile(optimizer=optimizer, loss='mean_squared_error', metrics=['accuracy'])
        return classifier

    nn_model_ps = KerasClassifier(build_fn=build_classifier)

    # Hyperparameters for the neural network
    para_nn = {
        'batch_size': [10, 32, 64],
        'epochs': [10, 100, 1000],
        'optimizer': ['adam', 'rmsprop', 'SGD'],
        'kernel': ['random_normal', 'random_uniform', 'truncated_normal'],
        'units': [8, 32, 64, 128],
        'hidden_layers': [2, 3, 5, 7],
        'activation': ['tanh', 'sigmoid', 'relu', 'selu']
    }

    # Randomized search for neural network tuning
    random_search_nn_ps = RandomizedSearchCV(estimator=nn_model_ps,
                                             param_distributions=para_nn,
                                             cv=skf.split(X_train, y_train),
                                             return_train_score=True,
                                             scoring='neg_brier_score', 
                                             n_iter=100)

    random_search_nn_ps.fit(cov_cols, t)

    # ------------------ XGBoost ------------------ #
    xgb_m = xgb.XGBClassifier(n_jobs=1, objective='binary:logistic', eval_metric="rmse")

    # Hyperparameters for XGBoost
    para_xgb = {
        'n_estimators': [1000, 600, 300, 100],
        'min_child_weight': [1, 10, 50],
        'gamma': [0.5, 2, None],
        'subsample': [0.6, 0.8, 1.0],
        'learning_rate': [0.02, 0.1, 0.2, 0.5],
        'max_depth': [3, 5, 7, 12]
    }

    # Randomized search for XGBoost hyperparameter tuning
    random_search_xgb_ps = RandomizedSearchCV(xgb_m, param_distributions=para_xgb,
                                              scoring='neg_brier_score',
                                              cv=skf.split(X_train, y_train),
                                              verbose=3, n_iter=100)

    random_search_xgb_ps.fit(cov_cols, t)

    # Predict propensity scores using all models
    PS_logcon = logmod_con.predict_proba(conf_cols)[:, 1]
    PS_lasso = logmod_ps.predict_proba(cov_cols)[:, 1]
    PS_nn_model = random_search_nn_ps.predict_proba(cov_cols)[:, 1]
    PS_xgb_m = random_search_xgb_ps.predict_proba(cov_cols)[:, 1]

    # Add propensity scores to covariates
    cov_cols['ps_log'] = PS_logcon
    cov_cols['ps_lasso'] = PS_lasso
    cov_cols['ps_xgb'] = PS_xgb_m
    cov_cols['ps_nn'] = PS_nn_model

    # Add treatment and outcome columns
    cov_cols['treatment'] = t
    cov_cols['Y'] = data["EVENT{}".format(i)]

    # Save the resulting data
    cov_cols.to_csv(output_file.format(i), index=False)


# Example of running the function
run_ps(
    i=1,
    conf_list=pd.read_csv("conf_list.csv", index_col=0),
    cov_list=pd.read_csv("cov_list.csv", index_col=0),
    data_template="D:/plasmode/data_{}.csv",
    output_file="D:/plasmode/output_m{}.csv"
)
