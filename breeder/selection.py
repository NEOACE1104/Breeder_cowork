"""Genomic prediction models implemented with pure Python linear algebra."""

from __future__ import annotations

import random

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .data import GenotypeDataset, PhenotypeDataset
from .utils import (
    k_fold_indices,
    matrix_inverse,
    matrix_multiply,
    matrix_vector_multiply,
    mean,
    pearson_correlation,
)


def _ridge_solution(x_matrix: List[List[float]], y_vector: List[float], alpha: float) -> List[float]:
    xt = list(zip(*x_matrix))
    xtx = matrix_multiply(xt, x_matrix)
    regularised = [[xtx[i][j] + (alpha if i == j else 0.0) for j in range(len(xtx))] for i in range(len(xtx))]
    xty = matrix_vector_multiply(xt, y_vector)
    inv = matrix_inverse(regularised)
    return matrix_vector_multiply(inv, xty)


def _standardise_matrix(matrix: List[List[float]]) -> tuple[List[List[float]], List[float], List[float]]:
    if not matrix:
        return [], [], []
    n_features = len(matrix[0])
    columns = list(zip(*matrix))
    means = [mean(col) for col in columns]
    stds: List[float] = []
    standardised: List[List[float]] = []
    for j in range(n_features):
        column = columns[j]
        mu = means[j]
        var = sum((value - mu) ** 2 for value in column)
        std = (var / (len(column) - 1)) ** 0.5 if len(column) > 1 else 0.0
        stds.append(std if std != 0 else 1.0)
    for row in matrix:
        standardised.append([(value - means[j]) / stds[j] if stds[j] else 0.0 for j, value in enumerate(row)])
    return standardised, means, stds


def _apply_standardisation(row: List[float], means: List[float], stds: List[float]) -> List[float]:
    return [(value - means[j]) / stds[j] if stds[j] else 0.0 for j, value in enumerate(row)]


def _rmse(true: List[float], pred: List[float]) -> float:
    if not true:
        return float("nan")
    return (sum((a - b) ** 2 for a, b in zip(true, pred)) / len(true)) ** 0.5


@dataclass
class GenomicBLUP:
    alpha: float = 1e-3
    scale: bool = True
    n_splits: int = 5
    random_state: Optional[int] = None
    marker_order_: List[str] = field(init=False, default_factory=list)
    mean_: List[float] = field(init=False, default_factory=list)
    scale_: List[float] = field(init=False, default_factory=list)
    intercept_: float | None = field(init=False, default=None)
    marker_effects_: List[float] | None = field(init=False, default=None)

    def _prepare(self, genotype: GenotypeDataset, phenotype: PhenotypeDataset) -> tuple[GenotypeDataset, List[float]]:
        pheno_map = phenotype.to_dict()
        common = [ind for ind in genotype.individuals if ind in pheno_map]
        aligned_geno = genotype.align(common)
        aligned_pheno = phenotype.subset(common)
        return aligned_geno, aligned_pheno.values

    def fit(self, genotype: GenotypeDataset, phenotype: PhenotypeDataset) -> "GenomicBLUP":
        geno, values = self._prepare(genotype, phenotype)
        matrix = geno.to_matrix(impute=True)
        if self.scale:
            matrix, means, stds = _standardise_matrix(matrix)
        else:
            means = [0.0 for _ in geno.markers]
            stds = [1.0 for _ in geno.markers]
        beta = _ridge_solution(matrix, values, self.alpha)
        self.marker_order_ = list(geno.markers)
        self.mean_ = means
        self.scale_ = stds
        self.intercept_ = mean(values)
        self.marker_effects_ = beta
        return self

    def predict(self, genotype: GenotypeDataset) -> List[float]:
        if self.marker_effects_ is None or self.intercept_ is None:
            raise RuntimeError("Model must be fitted before prediction.")
        reordered = genotype.reorder_markers(self.marker_order_)
        matrix = reordered.to_matrix(impute=True)
        if self.scale:
            matrix = [_apply_standardisation(row, self.mean_, self.scale_) for row in matrix]
        predictions = []
        for row in matrix:
            score = sum(effect * value for effect, value in zip(self.marker_effects_, row))
            predictions.append(self.intercept_ + score)
        return predictions

    def cross_validate(
        self,
        genotype: GenotypeDataset,
        phenotype: PhenotypeDataset,
        *,
        n_splits: Optional[int] = None,
        random_state: Optional[int] = None,
    ) -> List[Dict[str, float]]:
        geno, values = self._prepare(genotype, phenotype)
        matrix = geno.to_matrix(impute=True)
        n_splits = n_splits or self.n_splits
        rng = random.Random(random_state if random_state is not None else self.random_state)
        folds = k_fold_indices(len(matrix), n_splits, rng=rng)
        results: List[Dict[str, float]] = []
        for fold_id, (train_idx, test_idx) in enumerate(folds, start=1):
            x_train = [matrix[i] for i in train_idx]
            y_train = [values[i] for i in train_idx]
            x_test = [matrix[i] for i in test_idx]
            y_test = [values[i] for i in test_idx]
            if self.scale:
                x_train_std, means, stds = _standardise_matrix(x_train)
                x_test_std = [_apply_standardisation(row, means, stds) for row in x_test]
            else:
                x_train_std = x_train
                x_test_std = x_test
            beta = _ridge_solution(x_train_std, y_train, self.alpha)
            intercept = mean(y_train)
            preds = [intercept + sum(effect * value for effect, value in zip(beta, row)) for row in x_test_std]
            fold_rmse = _rmse(y_test, preds)
            corr = pearson_correlation(y_test, preds)
            results.append({"fold": float(fold_id), "rmse": fold_rmse, "correlation": corr, "n_test": float(len(test_idx))})
        return results
