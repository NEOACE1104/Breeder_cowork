import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency

def load_data_from_db(genotypes_df: pd.DataFrame, phenotypes_df: pd.DataFrame):
    """Pivots database data into wide format matrices."""
    geno_matrix = genotypes_df.pivot(index='sample_id', columns='snp_id', values='gt')
    pheno_matrix = phenotypes_df.pivot(index='sample_id', columns='trait', values='value')
    return geno_matrix, pheno_matrix

def impute_missing(geno_matrix: pd.DataFrame):
    """Imputes missing values (NaN) with the column mean. Fills with 0 if a column is all NaN."""
    geno_matrix_imputed = geno_matrix.apply(lambda x: x.fillna(x.mean()), axis=0)
    # If any columns are still all NaN (because they were all missing to begin with), fill them with 0
    geno_matrix_imputed.fillna(0, inplace=True)
    return geno_matrix_imputed

def filter_samples(geno_matrix: pd.DataFrame, max_sample_missing: float = 0.2):
    """Filters samples based on missing rate."""
    missing_rates = geno_matrix.isnull().mean(axis=1)
    keep_samples = missing_rates[missing_rates <= max_sample_missing].index
    return geno_matrix.loc[keep_samples], (len(geno_matrix) - len(keep_samples))

def filter_snps_by_missing(geno_matrix: pd.DataFrame, max_snp_missing: float = 0.2):
    """Filters SNPs based on missing rate."""
    missing_rates = geno_matrix.isnull().mean(axis=0)
    keep_snps = missing_rates[missing_rates <= max_snp_missing].index
    return geno_matrix[keep_snps], (len(geno_matrix.columns) - len(keep_snps))

def calculate_maf(geno_matrix: pd.DataFrame):
    """Calculates Minor Allele Frequency for each SNP."""
    p = geno_matrix.mean(axis=0) / 2
    maf = np.minimum(p, 1 - p)
    return maf

def filter_snps_by_maf(geno_matrix: pd.DataFrame, min_maf: float = 0.05):
    """Filters SNPs based on MAF."""
    maf = calculate_maf(geno_matrix)
    keep_snps = maf[maf >= min_maf].index
    return geno_matrix[keep_snps], (len(geno_matrix.columns) - len(keep_snps))

def hwe_test(genotypes: pd.Series):
    """
    Performs a Chi-squared test for Hardy-Weinberg Equilibrium.
    Returns the p-value.
    """
    genotypes = genotypes.dropna()
    counts = genotypes.value_counts()
    n0 = counts.get(0.0, 0)
    n1 = counts.get(1.0, 0)
    n2 = counts.get(2.0, 0)

    N = n0 + n1 + n2
    if N == 0:
        return 1.0

    # allele frequencies
    p_freq = (2 * n0 + n1) / (2 * N)
    q_freq = 1 - p_freq

    if p_freq == 0 or q_freq == 0:
        return 1.0

    # expected genotype counts
    exp_n0 = (p_freq**2) * N
    exp_n1 = (2 * p_freq * q_freq) * N
    exp_n2 = (q_freq**2) * N

    observed = [n0, n1, n2]
    expected = [exp_n0, exp_n1, exp_n2]

    if sum(observed) == 0: return 1.0

    # Chi-squared test
    stat, p_val, _, _ = chi2_contingency([observed, expected])
    return p_val if not np.isnan(p_val) else 1.0

def filter_snps_by_hwe(geno_matrix: pd.DataFrame, min_hwe_p: float = 1e-6):
    """Filters SNPs based on HWE p-value."""
    p_values = geno_matrix.apply(hwe_test, axis=0)
    keep_snps = p_values[p_values >= min_hwe_p].index
    return geno_matrix[keep_snps], (len(geno_matrix.columns) - len(keep_snps))

def run_qc_pipeline(geno_matrix: pd.DataFrame, params: dict):
    """Runs the full QC pipeline in the correct order."""
    log = []
    initial_samples, initial_snps = geno_matrix.shape
    log.append(f"Initial data: {initial_samples} samples, {initial_snps} SNPs.")

    # Step 1: Mark special missing values (0.0 -> NaN) based on context
    geno_processed = geno_matrix.copy()
    for snp in geno_processed.columns:
        counts = geno_processed[snp].value_counts()
        has_zero = 0.0 in counts.index
        has_others = any(x in counts.index for x in [1.0, 2.0])
        if has_zero and not has_others:
            geno_processed[snp] = geno_processed[snp].replace(0.0, np.nan)
    log.append("Step 1: Identified special 0.0 cases as missing values.")

    # Step 2: Sample filtering
    geno_processed, removed = filter_samples(geno_processed, params['max_sample_missing'])
    log.append(f"Step 2: Removed {removed} samples (missing rate > {params['max_sample_missing']}).")

    # Step 3: SNP filtering by missing rate
    geno_processed, removed = filter_snps_by_missing(geno_processed, params['max_snp_missing'])
    log.append(f"Step 3: Removed {removed} SNPs (missing rate > {params['max_snp_missing']}).")

    # Step 4: Impute remaining missing values
    geno_processed = impute_missing(geno_processed)
    log.append("Step 4: Imputed remaining missing values by SNP mean (or 0).")

    # Step 5: MAF filtering
    maf = calculate_maf(geno_processed)
    geno_processed, removed = filter_snps_by_maf(geno_processed, params['min_maf'])
    log.append(f"Step 5: Removed {removed} SNPs (MAF < {params['min_maf']}).")

    # Step 6: HWE filtering
    geno_processed, removed = filter_snps_by_hwe(geno_processed, params['min_hwe_p'])
    log.append(f"Step 6: Removed {removed} SNPs (HWE p-value < {params['min_hwe_p']}).")

    final_samples, final_snps = geno_processed.shape
    log.append(f"QC finished. Final data: {final_samples} samples, {final_snps} SNPs.")

    summary = {
        "initial_samples": initial_samples,
        "initial_snps": initial_snps,
        "samples_after_qc": final_samples,
        "snps_after_qc": final_snps,
        "mean_maf": maf.mean() if not maf.empty else 0.0,
    }

    return geno_processed, summary, "\n".join(log)
