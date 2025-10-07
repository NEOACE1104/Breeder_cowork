"""Analytical tools for genotype and phenotype data."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional

from .data import GenotypeDataset, PhenotypeDataset
from .utils import (
    covariance_matrix,
    matrix_inverse,
    matrix_vector_multiply,
    student_t_two_tailed_pvalue,
    symmetric_eigendecomposition,
)


@dataclass
class PopulationStructureAnalyzer:
    n_components: int = 10
    scale: bool = True
    components_: List[List[float]] = field(init=False, default_factory=list)
    eigenvalues_: List[float] = field(init=False, default_factory=list)
    means_: List[float] = field(init=False, default_factory=list)
    stds_: List[float] = field(init=False, default_factory=list)
    marker_order_: List[str] = field(init=False, default_factory=list)

    def fit(self, dataset: GenotypeDataset) -> "PopulationStructureAnalyzer":
        matrix, means, stds = dataset.standardized_matrix(impute=True)
        if not self.scale:
            matrix = dataset.to_matrix(impute=True)
            means = [0.0 for _ in dataset.markers]
            stds = [1.0 for _ in dataset.markers]
        cov = covariance_matrix(matrix)
        eigenvalues, eigenvectors = symmetric_eigendecomposition(cov, self.n_components)
        self.components_ = eigenvectors
        self.eigenvalues_ = eigenvalues
        self.means_ = means
        self.stds_ = stds
        self.marker_order_ = list(dataset.markers)
        return self

    def transform(self, dataset: GenotypeDataset) -> List[List[float]]:
        if not self.components_:
            raise RuntimeError("PopulationStructureAnalyzer must be fitted before transform().")
        reordered = dataset.reorder_markers(self.marker_order_)
        matrix: List[List[float]] = []
        for row in reordered.matrix:
            standardized_row = []
            for j, value in enumerate(row):
                val = value if value is not None else self.means_[j]
                standardized_row.append((val - self.means_[j]) / self.stds_[j] if self.stds_[j] else 0.0)
            matrix.append(standardized_row)
        scores: List[List[float]] = []
        for row in matrix:
            row_scores = []
            for component in self.components_:
                score = sum(row[j] * component[j] for j in range(len(component)))
                row_scores.append(score)
            scores.append(row_scores)
        return scores

    @property
    def explained_variance_ratio(self) -> List[float]:
        if not self.eigenvalues_:
            raise RuntimeError("Analyzer has not been fitted.")
        total = sum(self.eigenvalues_)
        if total == 0:
            return [0.0 for _ in self.eigenvalues_]
        return [value / total for value in self.eigenvalues_]


@dataclass
class AssociationAnalyzer:
    min_call_rate: float = 0.9
    impute: bool = True

    def run(
        self,
        genotype: GenotypeDataset,
        phenotype: PhenotypeDataset,
        covariates: Optional[Dict[str, List[float]]] = None,
    ) -> List[Dict[str, float]]:
        pheno_map = phenotype.to_dict()
        common = [ind for ind in genotype.individuals if ind in pheno_map]
        geno = genotype.align(common)
        pheno = phenotype.subset(common)
        cov_rows: List[List[float]] = []
        if covariates:
            cov_rows = [covariates.get(ind, []) for ind in geno.individuals]
        else:
            cov_rows = [[ ] for _ in geno.individuals]

        call_rates = geno.call_rate()
        results: List[Dict[str, float]] = []
        column_means, _ = geno.column_statistics()
        for marker_index, marker in enumerate(geno.markers):
            if call_rates[marker] < self.min_call_rate:
                continue
            design: List[List[float]] = []
            response: List[float] = []
            for i, individual in enumerate(geno.individuals):
                value = geno.matrix[i][marker_index]
                if value is None and not self.impute:
                    continue
                geno_value = value if value is not None else column_means[marker_index]
                row = [1.0]
                row.extend(cov_rows[i])
                row.append(geno_value)
                design.append(row)
                response.append(pheno.values[i])
            if len(design) <= len(design[0]):
                raise ValueError("Insufficient samples relative to covariates for marker analysis")
            Xt = list(zip(*design))
            XtX = [[sum(Xt[i][k] * design[k][j] for k in range(len(design))) for j in range(len(design[0]))] for i in range(len(Xt))]
            Xty = [sum(Xt[i][k] * response[k] for k in range(len(response))) for i in range(len(Xt))]
            XtX_inv = matrix_inverse(XtX)
            beta = matrix_vector_multiply(XtX_inv, Xty)
            predictions = [sum(beta[j] * design[i][j] for j in range(len(beta))) for i in range(len(design))]
            residuals = [response[i] - predictions[i] for i in range(len(response))]
            rss = sum(res ** 2 for res in residuals)
            dof = len(response) - len(beta)
            if dof <= 0:
                raise ValueError("Degrees of freedom must be positive")
            sigma2 = rss / dof
            variances = [XtX_inv[j][j] * sigma2 for j in range(len(beta))]
            effect = beta[-1]
            se = variances[-1] ** 0.5
            t_stat = effect / se if se > 0 else 0.0
            pvalue = student_t_two_tailed_pvalue(t_stat, dof)
            results.append(
                {
                    "marker": marker,
                    "effect": effect,
                    "se": se,
                    "t": t_stat,
                    "pvalue": pvalue,
                    "n": float(len(response)),
                }
            )
        return results
