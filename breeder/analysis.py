"""Analytical tools for genotype and phenotype data."""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .data import GenotypeDataset, PhenotypeDataset
from .utils import (
    covariance_matrix,
    matrix_inverse,
    matrix_multiply,
    transpose,
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


@dataclass
class DendrogramNode:
    """Node describing a hierarchical clustering result."""

    name: str
    distance: float = 0.0
    left: Optional["DendrogramNode"] = None
    right: Optional["DendrogramNode"] = None

    def is_leaf(self) -> bool:
        return self.left is None and self.right is None

    def members(self) -> List[str]:
        if self.is_leaf():
            return [self.name]
        members: List[str] = []
        if self.left is not None:
            members.extend(self.left.members())
        if self.right is not None:
            members.extend(self.right.members())
        return members


@dataclass
class DendrogramAnalyzer:
    """Perform simple agglomerative clustering on genotype data."""

    linkage: str = "average"
    root_: Optional[DendrogramNode] = field(init=False, default=None)
    merges_: List[Tuple[List[str], List[str], float]] = field(init=False, default_factory=list)

    def fit(self, dataset: GenotypeDataset) -> "DendrogramAnalyzer":
        matrix, _, _ = dataset.standardized_matrix(impute=True)
        n = len(dataset.individuals)
        if n == 0:
            raise ValueError("Dataset must contain at least one individual")

        clusters: Dict[int, DendrogramNode] = {
            i: DendrogramNode(name=dataset.individuals[i]) for i in range(n)
        }
        sizes: Dict[int, int] = {i: 1 for i in range(n)}
        distances: Dict[int, Dict[int, float]] = {i: {} for i in range(n)}

        def _euclidean(a: Sequence[float], b: Sequence[float]) -> float:
            return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))

        for i in range(n):
            for j in range(i + 1, n):
                dist = _euclidean(matrix[i], matrix[j])
                distances[i][j] = dist
                distances[j][i] = dist

        active_ids = list(clusters.keys())
        next_id = n
        self.merges_.clear()

        while len(active_ids) > 1:
            min_dist = float("inf")
            pair: Tuple[int, int] | None = None
            for i in active_ids:
                for j in active_ids:
                    if i >= j:
                        continue
                    dist = distances[i].get(j)
                    if dist is None:
                        continue
                    if dist < min_dist:
                        min_dist = dist
                        pair = (i, j)
            if pair is None:
                break
            i, j = pair
            left = clusters[i]
            right = clusters[j]
            new_node = DendrogramNode(
                name=f"cluster_{next_id}",
                distance=min_dist,
                left=left,
                right=right,
            )
            clusters[next_id] = new_node
            sizes[next_id] = sizes[i] + sizes[j]
            self.merges_.append((left.members(), right.members(), min_dist))

            for other in active_ids:
                if other in (i, j):
                    continue
                dist_i = distances[i].get(other, 0.0)
                dist_j = distances[j].get(other, 0.0)
                if self.linkage == "single":
                    new_dist = min(dist_i, dist_j)
                elif self.linkage == "complete":
                    new_dist = max(dist_i, dist_j)
                else:  # average linkage
                    new_dist = (
                        dist_i * sizes[i] + dist_j * sizes[j]
                    ) / (sizes[i] + sizes[j])
                distances.setdefault(other, {})[next_id] = new_dist
                distances.setdefault(next_id, {})[other] = new_dist

            for key in (i, j):
                active_ids.remove(key)
                distances.pop(key, None)
                sizes.pop(key, None)
                for dist_map in distances.values():
                    dist_map.pop(key, None)
                clusters.pop(key, None)

            active_ids.append(next_id)
            next_id += 1

        if not clusters:
            raise RuntimeError("Clustering failed to produce any nodes")
        self.root_ = clusters[active_ids[0]]
        return self

    def get_root(self) -> DendrogramNode:
        if self.root_ is None:
            raise RuntimeError("Analyzer has not been fitted")
        return self.root_


@dataclass
class AdmixtureAnalyzer:
    """Estimate ancestry proportions using a basic NMF approach."""

    n_populations: int
    max_iter: int = 500
    tol: float = 1e-4
    random_state: Optional[int] = None
    q_matrix_: List[List[float]] = field(init=False, default_factory=list)
    allele_frequencies_: List[List[float]] = field(init=False, default_factory=list)

    def fit(self, dataset: GenotypeDataset) -> "AdmixtureAnalyzer":
        if self.n_populations < 1:
            raise ValueError("n_populations must be at least 1")
        matrix = dataset.to_matrix(impute=True)
        if not matrix:
            raise ValueError("Genotype matrix is empty")
        n_individuals = len(matrix)
        n_markers = len(matrix[0])
        rng = random.Random(self.random_state)

        q = [[rng.random() + 1e-3 for _ in range(self.n_populations)] for _ in range(n_individuals)]
        h = [[rng.random() + 1e-3 for _ in range(n_markers)] for _ in range(self.n_populations)]

        def _normalize_rows(mat: List[List[float]]) -> None:
            for row in mat:
                total = sum(row)
                if total <= 0:
                    continue
                for idx in range(len(row)):
                    row[idx] = row[idx] / total

        _normalize_rows(q)

        last_error = float("inf")
        for _ in range(self.max_iter):
            # Update allele frequencies (H matrix)
            q_t = transpose(q)
            numerator_h = matrix_multiply(q_t, matrix)
            qtq = matrix_multiply(q_t, q)
            denominator_h = matrix_multiply(qtq, h)
            for k in range(self.n_populations):
                for j in range(n_markers):
                    denom = denominator_h[k][j] if denominator_h[k][j] > 1e-9 else 1e-9
                    h[k][j] = max(1e-9, h[k][j] * numerator_h[k][j] / denom)

            # Update q matrix (W matrix)
            h_t = transpose(h)
            numerator_q = matrix_multiply(matrix, h_t)
            hh_t = matrix_multiply(h, h_t)
            denominator_q = matrix_multiply(q, hh_t)
            for i in range(n_individuals):
                for k in range(self.n_populations):
                    denom = denominator_q[i][k] if denominator_q[i][k] > 1e-9 else 1e-9
                    q[i][k] = max(1e-9, q[i][k] * numerator_q[i][k] / denom)

            _normalize_rows(q)

            reconstruction = matrix_multiply(q, h)
            error = 0.0
            for i in range(n_individuals):
                for j in range(n_markers):
                    diff = matrix[i][j] - reconstruction[i][j]
                    error += diff * diff
            if abs(last_error - error) < self.tol:
                break
            last_error = error

        self.q_matrix_ = q
        self.allele_frequencies_ = h
        return self

    def transform(self, dataset: GenotypeDataset) -> List[List[float]]:
        if not self.q_matrix_:
            raise RuntimeError("AdmixtureAnalyzer must be fitted before transform().")
        # For simplicity, recompute memberships using fitted allele frequencies.
        matrix = dataset.to_matrix(impute=True)
        h_t = transpose(self.allele_frequencies_)
        hh_t = matrix_multiply(self.allele_frequencies_, h_t)
        q = [[1.0 for _ in range(self.n_populations)] for _ in matrix]
        numerator_q = matrix_multiply(matrix, h_t)
        denominator_q = matrix_multiply(q, hh_t)
        for i in range(len(matrix)):
            for k in range(self.n_populations):
                denom = denominator_q[i][k] if denominator_q[i][k] > 1e-9 else 1e-9
                q[i][k] = max(1e-9, q[i][k] * numerator_q[i][k] / denom)
            total = sum(q[i])
            if total > 0:
                q[i] = [value / total for value in q[i]]
        return q
