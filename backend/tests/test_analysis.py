import pytest
import pandas as pd
import numpy as np

from backend.analysis import qc, gwas

# --- Sample Data ---
def create_sample_data():
    """Creates a small, wide-format DataFrame for testing."""
    data = {
        'sample_id': [f's{i}' for i in range(10)],
        'trait_yield': np.random.rand(10) * 10,
        'snp_1_100': [0, 1, 2, 0, 1, 2, 0, 1, np.nan, 0], # Has one nan
        'snp_1_200': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], # Should be treated as all missing
        'snp_2_300': [0, 1, 1, 2, 2, 0, 0, 1, 1, 2], # Normal SNP
        'snp_2_400': [0, 0, 0, 0, 1, 0, 0, 0, 0, 0], # Low MAF (MAF = 0.05)
    }
    df = pd.DataFrame(data)

    pheno_cols = ['trait_yield']
    snp_cols = ['snp_1_100', 'snp_1_200', 'snp_2_300', 'snp_2_400']

    geno_matrix = df.set_index('sample_id')[snp_cols]
    pheno_matrix = df.set_index('sample_id')[pheno_cols]

    return geno_matrix, pheno_matrix

# --- QC Tests ---

def test_impute_missing():
    geno, _ = create_sample_data()
    # Manually set the "all zero" column to NaN to simulate first step of pipeline
    geno['snp_1_200'] = np.nan

    imputed_geno = qc.impute_missing(geno)

    # Check that no NaNs remain
    assert imputed_geno.isnull().sum().sum() == 0
    # Check that the all-NaN column was imputed to 0
    assert imputed_geno['snp_1_200'].mean() == 0

def test_filter_snps_by_maf():
    geno, _ = create_sample_data()
    # Impute NaNs first for a clean MAF calculation
    geno_imputed = qc.impute_missing(geno)

    filtered_geno, _ = qc.filter_snps_by_maf(geno_imputed, min_maf=0.1)

    assert 'snp_2_400' not in filtered_geno.columns
    assert 'snp_1_100' in filtered_geno.columns
    assert 'snp_2_300' in filtered_geno.columns

def test_qc_pipeline():
    geno, _ = create_sample_data()
    params = {
        'max_sample_missing': 0.3, # s8 has 2/4=0.5 missing. Others have 1/4=0.25. This keeps 9 samples.
        'max_snp_missing': 0.5,  # snp_1_200 is 100% missing, will be removed
        'min_maf': 0.1,          # snp_2_400 has MAF ~0.05, will be removed
        'min_hwe_p': 1e-6
    }

    qc_geno, summary, log = qc.run_qc_pipeline(geno, params)

    # s8 should be removed
    assert 's8' not in qc_geno.index
    assert qc_geno.shape[0] < geno.shape[0]

    # snp_1_200 and snp_2_400 should be removed
    assert 'snp_1_200' not in qc_geno.columns
    assert 'snp_2_400' not in qc_geno.columns
    assert qc_geno.shape[1] < geno.shape[1]

    assert summary['initial_snps'] == 4
    assert summary['initial_samples'] == 10
    assert summary['samples_after_qc'] == 9
    assert summary['snps_after_qc'] == 2 # snp_1_100 and snp_2_300 remain

# --- GWAS Tests ---

def test_run_pca():
    geno, _ = create_sample_data()
    # Run imputation to get a clean matrix for PCA
    geno_imputed = qc.impute_missing(geno)
    pcs = gwas.run_pca(geno_imputed, n_components=2)

    assert pcs.shape == (10, 2)
    assert 'PC1' in pcs.columns
    assert not pcs.isnull().values.any()

def test_gwas_ols_pipeline():
    geno, pheno = create_sample_data()
    # Impute to get clean data
    geno_imputed = qc.impute_missing(geno)

    request = {
        'trait': 'trait_yield',
        'method': 'ols',
        'n_pcs': 2
    }

    results_df = gwas.run_gwas_pipeline(geno_imputed, pheno, request)

    assert not results_df.empty
    assert 'p_value' in results_df.columns
    assert 'p_bonferroni' in results_df.columns
    # snp_1_200 has zero variance and will be skipped by statsmodels, so we expect 3 results, not 4
    assert len(results_df) == 3
    assert not results_df['p_value'].isnull().any()
