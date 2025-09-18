import matplotlib
matplotlib.use('Agg') # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from io import BytesIO
import base64
from scipy import stats

def get_image_as_base64(fig):
    """Converts a matplotlib figure to a base64 encoded string."""
    buf = BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight')
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode('utf-8')

def plot_maf_distribution(maf_series: pd.Series):
    """Plots a histogram of Minor Allele Frequencies."""
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.histplot(maf_series, bins=30, kde=True, ax=ax)
    ax.set_title('MAF Distribution')
    ax.set_xlabel('Minor Allele Frequency')
    ax.set_ylabel('Number of SNPs')
    return fig

def plot_manhattan(gwas_results: pd.DataFrame, marker_info: pd.DataFrame, bonferroni_thresh: float, fdr_thresh: float):
    """Generates a Manhattan plot for GWAS results."""
    df = gwas_results.join(marker_info, on='snp_id')
    df['-log10p'] = -np.log10(df['p_value'])
    df['chr'] = df['chr'].astype('category')
    df['chr_int'] = df['chr'].cat.codes

    # Sort by chromosome and position
    df = df.sort_values(['chr_int', 'bp'])
    df.reset_index(inplace=True, drop=True)
    df['ind'] = df.index

    fig, ax = plt.subplots(figsize=(16, 7))

    colors = ['#2c3e50', '#95a5a6']
    for i, (name, group) in enumerate(df.groupby('chr_int')):
        group.plot(kind='scatter', x='ind', y='-log10p', color=colors[i % len(colors)], ax=ax, s=10)

    ax.axhline(y=-np.log10(bonferroni_thresh), color='r', linestyle='--', label=f'Bonferroni (p={bonferroni_thresh:.2e})')
    ax.axhline(y=-np.log10(fdr_thresh), color='b', linestyle='--', label=f'FDR (p={fdr_thresh:.2e})')

    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title('Manhattan Plot')

    # Set chromosome ticks
    chrom_df = df.groupby('chr_int')['ind'].median()
    ax.set_xticks(chrom_df.values)
    ax.set_xticklabels(marker_info['chr'].unique())
    ax.legend()

    return fig

def plot_qq(gwas_results: pd.DataFrame):
    """Generates a QQ plot for GWAS p-values."""
    p_values = gwas_results['p_value'].dropna()

    observed_p = -np.log10(np.sort(p_values))
    expected_p = -np.log10(np.linspace(1/len(p_values), 1, len(p_values)))

    # Calculate lambda GC
    chisq = stats.chi2.ppf(1 - p_values, 1)
    lambda_gc = np.median(chisq) / stats.chi2.ppf(0.5, 1)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(expected_p, observed_p, s=10)
    ax.plot([0, max(expected_p)], [0, max(expected_p)], color='red', linestyle='--')
    ax.set_xlabel('Expected -log10(p)')
    ax.set_ylabel('Observed -log10(p)')
    ax.set_title('Q-Q Plot')
    ax.text(0.05, 0.9, f'λGC = {lambda_gc:.3f}', transform=ax.transAxes)

    return fig

def plot_mating_tradeoff(plot_data: dict):
    """Plots the trade-off between predicted gain and diversity."""
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(x=plot_data['diversities'], y=plot_data['gains'], ax=ax)
    ax.set_xlabel('Diversity Score (1 - Coancestry)')
    ax.set_ylabel('Predicted Genetic Gain')
    ax.set_title('Mating Pairs: Genetic Gain vs. Diversity')
    return fig
