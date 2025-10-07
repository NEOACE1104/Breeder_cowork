"""Core data structures for genotype and phenotype information."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .utils import mean, standard_deviation


def _parse_value(value: str) -> Optional[float]:
    value = value.strip()
    if value == "" or value.lower() == "nan":
        return None
    return float(value)


@dataclass
class GenotypeDataset:
    individuals: List[str]
    markers: List[str]
    matrix: List[List[Optional[float]]]

    @classmethod
    def from_csv(
        cls,
        path: str | Path,
        *,
        index_col: Optional[str] = None,
        delimiter: str = ",",
    ) -> "GenotypeDataset":
        with open(path, "r", newline="") as handle:
            reader = csv.reader(handle, delimiter=delimiter)
            header = next(reader)
            if index_col is None:
                index_idx = 0
            else:
                if index_col not in header:
                    raise KeyError(f"Index column '{index_col}' not found in genotype file")
                index_idx = header.index(index_col)
            markers = [value for i, value in enumerate(header) if i != index_idx]
            individuals: List[str] = []
            matrix: List[List[Optional[float]]] = []
            for row in reader:
                individuals.append(row[index_idx])
                values = [_parse_value(value) for i, value in enumerate(row) if i != index_idx]
                matrix.append(values)
        return cls(individuals, markers, matrix)

    def _column(self, marker_index: int) -> List[Optional[float]]:
        return [row[marker_index] for row in self.matrix]

    def allele_frequencies(self) -> Dict[str, float]:
        frequencies: Dict[str, float] = {}
        for j, marker in enumerate(self.markers):
            column = self._column(j)
            calls = sum(2 for value in column if value is not None)
            if calls == 0:
                frequencies[marker] = float("nan")
            else:
                allele_sum = sum(value for value in column if value is not None)
                frequencies[marker] = allele_sum / calls
        return frequencies

    def minor_allele_frequencies(self) -> Dict[str, float]:
        frequencies = self.allele_frequencies()
        return {marker: min(freq, 1 - freq) if not _is_nan(freq) else float("nan") for marker, freq in frequencies.items()}

    def call_rate(self) -> Dict[str, float]:
        rates: Dict[str, float] = {}
        total = len(self.individuals)
        for j, marker in enumerate(self.markers):
            non_missing = sum(1 for value in self._column(j) if value is not None)
            rates[marker] = non_missing / total if total else float("nan")
        return rates

    def heterozygosity(self) -> Dict[str, float]:
        values: Dict[str, float] = {}
        for j, marker in enumerate(self.markers):
            column = self._column(j)
            het = sum(1 for value in column if value == 1)
            total = sum(1 for value in column if value is not None)
            values[marker] = het / total if total else float("nan")
        return values

    def column_statistics(self) -> Tuple[List[float], List[float]]:
        means: List[float] = []
        stds: List[float] = []
        for j in range(len(self.markers)):
            column = [value for value in self._column(j) if value is not None]
            if column:
                col_mean = mean(column)
                col_std = standard_deviation(column, ddof=1) if len(column) > 1 else 0.0
            else:
                col_mean = 0.0
                col_std = 0.0
            means.append(col_mean)
            stds.append(col_std if col_std != 0 else 1.0)
        return means, stds

    def to_matrix(self, *, impute: bool) -> List[List[float]]:
        means, _ = self.column_statistics()
        filled: List[List[float]] = []
        for row in self.matrix:
            filled.append([value if value is not None else means[j] for j, value in enumerate(row)] if impute else [float(value) if value is not None else float("nan") for value in row])
        return filled

    def standardized_matrix(self, *, impute: bool = True) -> Tuple[List[List[float]], List[float], List[float]]:
        means, stds = self.column_statistics()
        matrix = []
        for row in self.matrix:
            standardized_row: List[float] = []
            for j, value in enumerate(row):
                val = value if value is not None else (means[j] if impute else 0.0)
                standardized_row.append((val - means[j]) / stds[j] if stds[j] else 0.0)
            matrix.append(standardized_row)
        return matrix, means, stds

    def reorder_markers(self, markers: Sequence[str]) -> "GenotypeDataset":
        index = {marker: idx for idx, marker in enumerate(self.markers)}
        new_matrix: List[List[Optional[float]]] = []
        for row in self.matrix:
            new_matrix.append([row[index[m]] if m in index else None for m in markers])
        return GenotypeDataset(list(self.individuals), list(markers), new_matrix)

    def align(self, individuals: Iterable[str]) -> "GenotypeDataset":
        individual_set = set(individuals)
        new_individuals: List[str] = []
        new_matrix: List[List[Optional[float]]] = []
        for name, row in zip(self.individuals, self.matrix):
            if name in individual_set:
                new_individuals.append(name)
                new_matrix.append(list(row))
        return GenotypeDataset(new_individuals, list(self.markers), new_matrix)


@dataclass
class PhenotypeDataset:
    individuals: List[str]
    values: List[float]

    @classmethod
    def from_csv(
        cls,
        path: str | Path,
        *,
        trait: str,
        index_col: Optional[str] = None,
        delimiter: str = ",",
    ) -> "PhenotypeDataset":
        with open(path, "r", newline="") as handle:
            reader = csv.reader(handle, delimiter=delimiter)
            header = next(reader)
            if index_col is None:
                index_idx = 0
            else:
                if index_col not in header:
                    raise KeyError(f"Index column '{index_col}' not found in phenotype file")
                index_idx = header.index(index_col)
            if trait not in header:
                raise KeyError(f"Trait '{trait}' not found in phenotype file")
            trait_idx = header.index(trait)
            individuals: List[str] = []
            values: List[float] = []
            for row in reader:
                value = _parse_value(row[trait_idx])
                if value is None:
                    continue
                individuals.append(row[index_idx])
                values.append(float(value))
        return cls(individuals, values)

    def to_dict(self) -> Dict[str, float]:
        return dict(zip(self.individuals, self.values))

    def summary(self) -> Dict[str, float]:
        ordered = sorted(self.values)
        n = len(ordered)
        if n == 0:
            raise ValueError("Phenotype dataset is empty")
        mid = n // 2
        median = (ordered[mid - 1] + ordered[mid]) / 2 if n % 2 == 0 else ordered[mid]
        q1 = ordered[n // 4]
        q3 = ordered[(3 * n) // 4]
        return {
            "count": float(n),
            "mean": mean(self.values),
            "std": standard_deviation(self.values, ddof=1) if n > 1 else 0.0,
            "min": ordered[0],
            "25%": q1,
            "50%": median,
            "75%": q3,
            "max": ordered[-1],
        }

    def normalize(self) -> "PhenotypeDataset":
        if not self.values:
            raise ValueError("No phenotype values to normalise")
        mu = mean(self.values)
        sd = standard_deviation(self.values, ddof=1) if len(self.values) > 1 else 0.0
        if sd == 0:
            raise ValueError("Phenotype variance is zero; cannot normalize")
        normalized = [(value - mu) / sd for value in self.values]
        return PhenotypeDataset(list(self.individuals), normalized)

    def subset(self, individuals: Iterable[str]) -> "PhenotypeDataset":
        selected = list(individuals)
        index = {ind: idx for idx, ind in enumerate(self.individuals)}
        values = [self.values[index[name]] for name in selected if name in index]
        return PhenotypeDataset(selected[: len(values)], values)

    def align(self, genotype: GenotypeDataset) -> "PhenotypeDataset":
        return self.subset(genotype.individuals)


def _is_nan(value: float) -> bool:
    return value != value
