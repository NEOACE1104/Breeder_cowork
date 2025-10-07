"""Utility functions for statistical computations without external dependencies."""

from __future__ import annotations

import math
import random
from typing import List, Sequence, Tuple


def mean(values: Sequence[float]) -> float:
    if not values:
        raise ValueError("mean requires at least one value")
    return sum(values) / len(values)


def variance(values: Sequence[float], *, ddof: int = 1) -> float:
    if len(values) <= ddof:
        raise ValueError("variance requires more values than degrees of freedom")
    mu = mean(values)
    return sum((x - mu) ** 2 for x in values) / (len(values) - ddof)


def standard_deviation(values: Sequence[float], *, ddof: int = 1) -> float:
    return math.sqrt(variance(values, ddof=ddof))


def transpose(matrix: Sequence[Sequence[float]]) -> List[List[float]]:
    if not matrix:
        return []
    return [list(col) for col in zip(*matrix)]


def matrix_multiply(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]]) -> List[List[float]]:
    if not a or not b:
        return []
    num_rows = len(a)
    num_cols = len(b[0])
    num_inner = len(b)
    result = [[0.0 for _ in range(num_cols)] for _ in range(num_rows)]
    for i in range(num_rows):
        for k in range(num_inner):
            aik = a[i][k]
            if aik == 0:
                continue
            for j in range(num_cols):
                result[i][j] += aik * b[k][j]
    return result


def matrix_vector_multiply(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> List[float]:
    return [sum(row[j] * vector[j] for j in range(len(vector))) for row in matrix]


def vector_dot(a: Sequence[float], b: Sequence[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def vector_norm(vector: Sequence[float]) -> float:
    return math.sqrt(vector_dot(vector, vector))


def normalize_vector(vector: Sequence[float]) -> List[float]:
    norm = vector_norm(vector)
    if norm == 0:
        raise ValueError("Cannot normalize zero vector")
    return [x / norm for x in vector]


def outer_product(a: Sequence[float], b: Sequence[float]) -> List[List[float]]:
    return [[x * y for y in b] for x in a]


def subtract_matrices(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]]) -> List[List[float]]:
    return [[x - y for x, y in zip(row_a, row_b)] for row_a, row_b in zip(a, b)]


def identity_matrix(n: int) -> List[List[float]]:
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def gauss_jordan(matrix: Sequence[Sequence[float]], rhs: Sequence[Sequence[float]]) -> List[List[float]]:
    n = len(matrix)
    aug = [list(row) + list(right_row) for row, right_row in zip(matrix, rhs)]
    m_cols = len(aug[0])
    for col in range(n):
        pivot_row = max(range(col, n), key=lambda r: abs(aug[r][col]))
        pivot = aug[pivot_row][col]
        if abs(pivot) < 1e-12:
            raise ValueError("Matrix is singular")
        aug[col], aug[pivot_row] = aug[pivot_row], aug[col]
        pivot = aug[col][col]
        for j in range(col, m_cols):
            aug[col][j] /= pivot
        for row in range(n):
            if row == col:
                continue
            factor = aug[row][col]
            if factor == 0:
                continue
            for j in range(col, m_cols):
                aug[row][j] -= factor * aug[col][j]
    return [row[n:] for row in aug]


def matrix_inverse(matrix: Sequence[Sequence[float]]) -> List[List[float]]:
    identity = identity_matrix(len(matrix))
    return gauss_jordan(matrix, identity)


def solve_linear_system(matrix: Sequence[Sequence[float]], rhs: Sequence[float]) -> List[float]:
    solutions = gauss_jordan(matrix, [[value] for value in rhs])
    return [row[0] for row in solutions]


def covariance_matrix(matrix: Sequence[Sequence[float]]) -> List[List[float]]:
    if not matrix:
        return []
    n_samples = len(matrix)
    n_features = len(matrix[0])
    means = [mean([matrix[i][j] for i in range(n_samples)]) for j in range(n_features)]
    centered = [[matrix[i][j] - means[j] for j in range(n_features)] for i in range(n_samples)]
    scaled = 1.0 / (n_samples - 1)
    cov = [[0.0 for _ in range(n_features)] for _ in range(n_features)]
    for i in range(n_features):
        for j in range(i, n_features):
            value = scaled * sum(centered[k][i] * centered[k][j] for k in range(n_samples))
            cov[i][j] = value
            cov[j][i] = value
    return cov


def power_iteration(matrix: Sequence[Sequence[float]], *, iterations: int = 1000, tol: float = 1e-9) -> Tuple[float, List[float]]:
    n = len(matrix)
    vector = [1.0 / math.sqrt(n) for _ in range(n)]
    last_value = 0.0
    for _ in range(iterations):
        next_vector = matrix_vector_multiply(matrix, vector)
        norm = vector_norm(next_vector)
        if norm == 0:
            break
        next_vector = [x / norm for x in next_vector]
        eigenvalue = vector_dot(next_vector, matrix_vector_multiply(matrix, next_vector))
        if abs(eigenvalue - last_value) < tol:
            vector = next_vector
            last_value = eigenvalue
            break
        vector = next_vector
        last_value = eigenvalue
    return last_value, vector


def symmetric_eigendecomposition(matrix: Sequence[Sequence[float]], n_components: int) -> Tuple[List[float], List[List[float]]]:
    residual = [list(row) for row in matrix]
    eigenvalues: List[float] = []
    eigenvectors: List[List[float]] = []
    for _ in range(min(n_components, len(matrix))):
        eigenvalue, eigenvector = power_iteration(residual)
        eigenvalues.append(eigenvalue)
        eigenvectors.append(eigenvector)
        outer = outer_product(eigenvector, eigenvector)
        scaled = [[eigenvalue * outer[i][j] for j in range(len(matrix))] for i in range(len(matrix))]
        residual = subtract_matrices(residual, scaled)
    return eigenvalues, eigenvectors


def pearson_correlation(x: Sequence[float], y: Sequence[float]) -> float:
    if len(x) != len(y):
        raise ValueError("Sequences must be of equal length")
    if len(x) < 2:
        return float("nan")
    mean_x = mean(x)
    mean_y = mean(y)
    num = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y))
    den_x = math.sqrt(sum((a - mean_x) ** 2 for a in x))
    den_y = math.sqrt(sum((b - mean_y) ** 2 for b in y))
    if den_x == 0 or den_y == 0:
        return float("nan")
    return num / (den_x * den_y)


def student_t_pdf(x: float, df: int) -> float:
    num = math.gamma((df + 1) / 2)
    den = math.sqrt(df * math.pi) * math.gamma(df / 2)
    return (num / den) * (1 + (x**2) / df) ** (-(df + 1) / 2)


def student_t_cdf(x: float, df: int) -> float:
    if x == 0:
        return 0.5
    if x < 0:
        return 1 - student_t_cdf(-x, df)
    upper = max(x, 1.0)
    steps = max(200, int(upper * 200))
    h = x / steps
    total = student_t_pdf(0.0, df) + student_t_pdf(x, df)
    for i in range(1, steps, 2):
        total += 4 * student_t_pdf(i * h, df)
    for i in range(2, steps, 2):
        total += 2 * student_t_pdf(i * h, df)
    area = total * h / 3
    return 0.5 + area


def student_t_two_tailed_pvalue(t_stat: float, df: int) -> float:
    prob = student_t_cdf(abs(t_stat), df)
    return max(0.0, min(1.0, 2 * (1 - prob)))


def chi_square_sf_1df(chi2: float) -> float:
    return math.erfc(math.sqrt(chi2 / 2))


def k_fold_indices(n_samples: int, n_splits: int, *, rng: random.Random) -> List[Tuple[List[int], List[int]]]:
    indices = list(range(n_samples))
    rng.shuffle(indices)
    fold_sizes = [n_samples // n_splits] * n_splits
    for i in range(n_samples % n_splits):
        fold_sizes[i] += 1
    folds = []
    current = 0
    for fold_size in fold_sizes:
        test_indices = indices[current : current + fold_size]
        train_indices = [idx for idx in indices if idx not in test_indices]
        folds.append((train_indices, test_indices))
        current += fold_size
    return folds
