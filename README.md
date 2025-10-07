# Breeder_cowork

Breeder is a scientific analysis toolkit for genomic breeding experiments. It provides
robust data structures, quality-control routines, population structure analysis, association
mapping, and genomic prediction models suitable for research-scale studies.

## Features

- **Flexible data ingestion** for genotype dosage matrices and phenotype traits.
- **Quality control metrics** including call rate, minor allele frequency, heterozygosity, and
  Hardy–Weinberg equilibrium testing with configurable thresholds.
- **Population structure characterisation** using PCA with optional scaling.
- **Genome-wide association analysis** via ordinary least squares regression supporting
  covariates, missing-data handling, and per-marker statistics.
- **Genomic BLUP predictions** with ridge regularisation, reusable marker effects, and
  cross-validation diagnostics.
- **Command line interface** for reproducible pipelines (QC, PCA, GWAS, GBLUP).

## Installation

The project is packaged as a Python module. Install the dependencies with:

```bash
pip install -e .
```

For development, include the optional test dependencies:

```bash
pip install -e .[tests]
```

## Command line usage

All tools are exposed through the `breeder` CLI.

### Quality control

```bash
breeder qc \
  --genotype path/to/genotypes.csv \
  --geno-index id \
  --metrics-out qc_metrics.csv \
  --filtered-genotype filtered_genotypes.csv
```

### Principal component analysis

```bash
breeder pca \
  --genotype path/to/genotypes.csv \
  --geno-index id \
  --components 5 \
  --output pcs.csv \
  --variance-out pca_variance.csv
```

### Genome-wide association study

```bash
breeder gwas \
  --genotype path/to/genotypes.csv \
  --geno-index id \
  --phenotype path/to/phenotypes.csv \
  --trait yield \
  --pheno-index id \
  --covariates path/to/covariates.csv \
  --cov-index id \
  --output gwas_results.csv
```

### Genomic BLUP

Cross-validate model performance:

```bash
breeder gblup \
  --genotype path/to/genotypes.csv \
  --geno-index id \
  --phenotype path/to/phenotypes.csv \
  --trait yield \
  --pheno-index id \
  --cross-validate \
  --output gblup_cv.csv
```

Fit a model and save predictions:

```bash
breeder gblup \
  --genotype path/to/genotypes.csv \
  --geno-index id \
  --phenotype path/to/phenotypes.csv \
  --trait yield \
  --pheno-index id \
  --output training_predictions.csv \
  --predict-genotype path/to/validation_genotypes.csv \
  --predict-index id \
  --predict-output validation_predictions.csv
```

## Testing

Run the unit test suite with:

```bash
pytest
```
