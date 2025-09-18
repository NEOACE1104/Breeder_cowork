import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.preprocessing import StandardScaler

def calculate_grm(X: pd.DataFrame):
    """Calculates the Genomic Relationship Matrix (GRM)."""
    # Standardize X: (X - mu) / sd
    scaler = StandardScaler()
    X_std = scaler.fit_transform(X)

    m = X_std.shape[1]
    G = (X_std @ X_std.T) / m
    return pd.DataFrame(G, index=X.index, columns=X.index)

def train_gblup_model(X_train: pd.DataFrame, y_train: pd.Series, n_splits=5):
    """
    Trains a GBLUP/rrBLUP model and evaluates it using K-fold cross-validation.
    """
    # For rrBLUP on SNP effects, we use Ridge regression on the individuals, which is equivalent to GBLUP.
    # The alpha parameter for Ridge corresponds to the regularization parameter lambda.
    # A common heuristic for lambda is `n_snps / variance_ratio`, but we can let RidgeCV find a good one.
    # For simplicity in MVP, we'll use a fixed alpha or a simple RidgeCV.

    model = Ridge(alpha=1.0) # A default alpha, can be tuned.

    # Cross-validation
    cv = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    y_pred_cv = cross_val_predict(model, X_train, y_train, cv=cv)

    # Calculate metrics
    r2 = r2_score(y_train, y_pred_cv)
    rmse = np.sqrt(mean_squared_error(y_train, y_pred_cv))
    pearson_r = np.corrcoef(y_train, y_pred_cv)[0, 1]

    cv_results = {
        "r2": r2,
        "rmse": rmse,
        "pearson_r": pearson_r,
    }

    # Fit final model on all training data
    model.fit(X_train, y_train)

    return model, cv_results

def predict_gebv(model: Ridge, X_predict: pd.DataFrame):
    """Predicts GEBVs for a new set of individuals."""
    gebvs = model.predict(X_predict)
    return pd.Series(gebvs, index=X_predict.index)

def run_gp_pipeline(geno_matrix: pd.DataFrame, pheno_matrix: pd.DataFrame, request: dict):
    """
    Full pipeline for a single trait GP.
    """
    trait = request['traits'][0] # MVP: single trait first
    train_ids = request['train_ids']
    candidate_ids = request.get('candidate_ids', [])

    # Align data
    common_samples = geno_matrix.index.intersection(pheno_matrix.index).intersection(train_ids)
    X_train = geno_matrix.loc[common_samples]
    y_train = pheno_matrix.loc[common_samples, trait]

    if X_train.empty or y_train.empty:
        raise ValueError("No valid training data found for the selected samples and trait.")

    # Train model
    model, cv_results = train_gblup_model(X_train, y_train, n_splits=request.get('k_folds', 5))

    # Predict candidates if any
    gebvs = None
    if candidate_ids:
        predict_ids_in_geno = list(set(candidate_ids) & set(geno_matrix.index))
        if predict_ids_in_geno:
            X_predict = geno_matrix.loc[predict_ids_in_geno]
            gebv_series = predict_gebv(model, X_predict)
            gebvs = gebv_series.reset_index().rename(columns={'index': 'sample_id', 0: 'gebv'}).to_dict('records')

    return model, cv_results, gebvs
