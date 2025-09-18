import pandas as pd
import numpy as np
import statsmodels.api as sm
from sklearn.decomposition import PCA

# Placeholder for limix, if available
try:
    from limix.qc import quantile_gaussianize
    from limix.qtl import qtl_test_lmm
    LIMIX_AVAILABLE = True
except ImportError:
    LIMIX_AVAILABLE = False

def run_pca(X: pd.DataFrame, n_components: int):
    """Performs PCA and returns the top components."""
    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(X)
    return pd.DataFrame(pcs, index=X.index, columns=[f'PC{i+1}' for i in range(n_components)])

def run_gwas_ols(
    geno_matrix: pd.DataFrame,
    pheno_series: pd.Series,
    covariates: pd.DataFrame
):
    """
    Performs GWAS using an Ordinary Least Squares (OLS) model for each SNP.
    y ~ snp + covariates
    """
    results = []

    # Align data
    common_samples = geno_matrix.index.intersection(pheno_series.index).intersection(covariates.index)
    geno_matrix = geno_matrix.loc[common_samples]
    pheno_series = pheno_series.loc[common_samples]
    covariates = covariates.loc[common_samples]

    X_cov = sm.add_constant(covariates)

    for snp_id in geno_matrix.columns:
        snp_vector = geno_matrix[snp_id]
        X_full = X_cov.copy()
        X_full['snp'] = snp_vector

        try:
            model = sm.OLS(pheno_series, X_full).fit()
            p_value = model.pvalues.get('snp', np.nan)
            effect_size = model.params.get('snp', np.nan)

            if not np.isnan(p_value):
                results.append({
                    'snp_id': snp_id,
                    'p_value': p_value,
                    'effect_size': effect_size
                })
        except Exception as e:
            # Handle cases where model fitting fails (e.g., perfect collinearity)
            print(f"Could not fit model for SNP {snp_id}: {e}")
            continue

    return pd.DataFrame(results)

def run_gwas_lmm(
    geno_matrix: pd.DataFrame,
    pheno_series: pd.Series,
    covariates: pd.DataFrame,
    kinship_matrix: pd.DataFrame
):
    """
    Performs GWAS using a Linear Mixed Model (LMM).
    This is a placeholder and requires limix.
    """
    if not LIMIX_AVAILABLE:
        raise ImportError("limix is not installed. LMM is not available.")

    # Align data
    common_samples = geno_matrix.index.intersection(pheno_series.index).intersection(covariates.index).intersection(kinship_matrix.index)
    geno_matrix = geno_matrix.loc[common_samples]
    pheno_series = pheno_series.loc[common_samples]
    covariates = covariates.loc[common_samples]
    kinship_matrix = kinship_matrix.loc[common_samples, common_samples]

    # Quantile normalize the phenotype
    pheno_series = quantile_gaussianize(pheno_series.values)

    # Run LMM
    # The `qtl_test_lmm` expects genotypes (samples x variants), covariates, and kinship matrix K
    lmm_results = qtl_test_lmm(geno_matrix.values, pheno_series, K=kinship_matrix.values, M=covariates.values, verbose=False)

    results_df = pd.DataFrame({
        'snp_id': geno_matrix.columns,
        'p_value': lmm_results.p_value,
        'effect_size': lmm_results.beta,
    })

    return results_df


def run_gwas_pipeline(geno_matrix: pd.DataFrame, pheno_matrix: pd.DataFrame, request: dict):
    """
    Full pipeline for GWAS analysis.
    """
    trait = request['trait']
    method = request.get('method', 'ols')
    n_pcs = request.get('n_pcs', 5)

    if trait not in pheno_matrix.columns:
        raise ValueError(f"Trait '{trait}' not found in phenotype data.")

    pheno_series = pheno_matrix[trait].dropna()

    # Align geno and pheno before PCA
    common_samples = geno_matrix.index.intersection(pheno_series.index)
    geno_aligned = geno_matrix.loc[common_samples]
    pheno_aligned = pheno_series.loc[common_samples]

    # Calculate PCs for population structure
    if n_pcs > 0:
        covariates = run_pca(geno_aligned, n_pcs)
    else:
        covariates = pd.DataFrame(index=geno_aligned.index)

    if method == 'lmm' and LIMIX_AVAILABLE:
        # from backend.analysis.gp import calculate_grm
        # kinship = calculate_grm(geno_aligned)
        # results_df = run_gwas_lmm(geno_aligned, pheno_aligned, covariates, kinship)
        print("LMM not fully implemented in this MVP, falling back to OLS.")
        results_df = run_gwas_ols(geno_aligned, pheno_aligned, covariates)
    else:
        results_df = run_gwas_ols(geno_aligned, pheno_aligned, covariates)

    # Add Bonferroni correction
    num_tests = len(results_df)
    results_df['p_bonferroni'] = np.minimum(1.0, results_df['p_value'] * num_tests)

    # Add FDR correction
    results_df = results_df.sort_values(by='p_value')
    results_df['p_fdr'] = results_df['p_value'] * num_tests / np.arange(1, num_tests + 1)
    results_df['p_fdr'] = np.minimum(1.0, results_df['p_fdr'].cummin()) # Ensure monotonicity

    return results_df.sort_index()
