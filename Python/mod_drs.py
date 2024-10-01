import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from keras.models import Sequential
from keras.layers import Dense
from keras.wrappers.scikit_learn import KerasClassifier
from sklearn.model_selection import StratifiedKFold
import xgboost as xgb
import copy


# Refactor into a function with customizable inputs
def run_drs(i, conf_list_path, cov_list_path, data_path_template, output_path_template, random_seed=43):
      """
    Function to run disease risk scores scores simulation on covariates, confounders, and treatment models
    
    Parameters:
    i (int): Simulation index
    conf_list_path (str): Path to confounder list
    cov_list_path (str): Path to covariate list
    data_path_template (str): Path template to load the data for each iteration
    output_path_template (str): Path template to save the output CSV
    random_seed (int): Random seed for reproducibility
    """
    # Load confounder and covariate lists
    conf_list = pd.read_csv(conf_list_path, index_col=0)
    cov_list = pd.read_csv(cov_list_path, index_col=0)

    # Load dataset for the current simulation based on the iteration index (i)
    simdata = pd.read_csv(data_path_template.format(i), index_col=0).dropna()

    # Create a deep copy of the simulation data
    data = pd.DataFrame(copy.deepcopy(simdata))

    # Filter confounders and covariates based on provided lists
    conf_cols = data.filter(map(str, conf_list['x'].to_list()))
    conf_cols['treatment'] = data["EXPOSURE{}".format(i)]

    cov_cols = data.filter(map(str, cov_list['x'].to_list()))
    cov_cols['treatment'] = data["EXPOSURE{}".format(i)]

    # Define treatment and outcome variables
    t = data["EXPOSURE{}".format(i)]
    Y = data["EVENT{}".format(i)]

    # Split data into training and testing sets (80% train, 20% test), stratified by treatment
    X_train, X_test, y_train, y_test = train_test_split(cov_cols, Y, test_size=0.8, stratify=t, random_state=random_seed)

    # Logistic regression model (confounders only)
    logmod_con = LogisticRegression(penalty='none', solver='lbfgs', max_iter=200)
    logmod_con.fit(conf_cols, Y)

    # Logistic regression with L1 regularization (LASSO)
    logistic = LogisticRegression(penalty='l1', max_iter=300)
    alphas = np.arange(0.01, 0.51, 0.05)
    hyperparameters = dict(C=alphas, solver=['liblinear', 'saga'])

    # Randomized search for hyperparameter tuning (LASSO)
    clf_ps = RandomizedSearchCV(estimator=logistic, param_distributions=hyperparameters, cv=5, verbose=0,
                                scoring='neg_brier_score', n_iter=15)
    logmod_ps = clf_ps.fit(X_train, y_train)

    # Fit the tuned LASSO model on full covariate data
    logmod_ps.fit(cov_cols, Y)

    # ------------------ Multi-Layer Perceptron (MLP) ------------------ #
    def build_classifier(optimizer, kernel, units, hidden_layers, activation):
        classifier = Sequential()
        classifier.add(Dense(units=units, activation='sigmoid', input_dim=cov_cols.shape[1], kernel_initializer=kernel))

        for _ in range(hidden_layers):
            classifier.add(Dense(units=units, activation=activation, kernel_initializer=kernel))

        classifier.add(Dense(1, activation='sigmoid', kernel_initializer=kernel))
        classifier.compile(optimizer=optimizer, loss='mean_squared_error', metrics=['accuracy'])
        return classifier

    nn_model_ps = KerasClassifier(build_fn=build_classifier)

    # Cross-validation and hyperparameter tuning setup for MLP
    para_nn = {'batch_size': [5, 8, 10],
               'epochs': [1, 10],
               'optimizer': ['adam', 'rmsprop', 'SGD'],
               'kernel': ['random_normal'],
               'units': [4, 8, 32],
               "hidden_layers": [2, 3, 5],
               "activation": ['tanh', 'sigmoid', 'relu', 'selu']}

    # Randomized search for neural network hyperparameter tuning
    random_search_nn_ps = RandomizedSearchCV(estimator=nn_model_ps, param_distributions=para_nn,
                                             cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=1001).split(X_train, y_train),
                                             return_train_score=True, scoring='neg_brier_score', n_iter=15)

    random_search_nn_ps.fit(cov_cols, Y)

    # ------------------ XGBoost ------------------ #
    xgb_m = xgb.XGBClassifier(n_jobs=1, objective='binary:logistic', eval_metric="mae")

    # Hyperparameters for XGBoost
    para_xgb = {
        'n_estimators': [600, 300, 100],
        'min_child_weight': [1, 10, 50],
        'gamma': [0.5, 2, None],
        'subsample': [0.6, 0.8, 1.0],
        'learning_rate': [0.02, 0.1, 0.2, 0.5],
        'max_depth': [3, 5, 7, 12]
    }

    # Randomized search for XGBoost hyperparameter tuning
    random_search_xgb_ps = RandomizedSearchCV(xgb_m, param_distributions=para_xgb, scoring='neg_brier_score',
                                              cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=1001).split(X_train, y_train),
                                              verbose=3, n_iter=15)

    random_search_xgb_ps.fit(cov_cols, Y)

    # Predict propensity scores using different models
    PS_logcon = logmod_con.predict_proba(conf_cols)[:, 1]
    PS_lasso = logmod_ps.predict_proba(cov_cols)[:, 1]
    PS_nn_model = random_search_nn_ps.predict_proba(cov_cols)[:, 1]
    PS_xgb_m = random_search_xgb_ps.predict_proba(cov_cols)[:, 1]

    # Store propensity scores in covariates
    cov_cols['ps_log'] = PS_logcon
    cov_cols['ps_lasso'] = PS_lasso
    cov_cols['ps_xgb'] = PS_xgb_m
    cov_cols['ps_nn'] = PS_nn_model

    # Add treatment and outcome columns to the covariates
    cov_cols['treatment'] = t
    cov_cols['Y'] = data["EVENT{}".format(i)]

    # Save the output data to CSV
    cov_cols.to_csv(output_path_template.format(i), index=False)


# Example of running the function
run_drs(
    i=30,
    conf_list_path="...",
    cov_list_path="...",
    data_path_template="...",
    output_path_template="..."
)
