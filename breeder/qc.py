"""Quality control routines for genotype datasets."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

from .data import GenotypeDataset
from .utils import chi_square_sf_1df


def hardy_weinberg_pvalue(column: List[Optional[float]]) -> float:
    counts = {0: 0, 1: 0, 2: 0}
    for value in column:
        if value is None:
            continue
        rounded = int(round(value))
        if rounded in counts:
            counts[rounded] += 1
    n = sum(counts.values())
    if n == 0:
        return float("nan")
    p = (2 * counts[2] + counts[1]) / (2 * n)
    q = 1 - p
    expected = [q**2 * n, 2 * p * q * n, p**2 * n]
    chi2 = 0.0
    for observed, exp in zip([counts[0], counts[1], counts[2]], expected):
        if exp == 0:
            continue
        chi2 += (observed - exp) ** 2 / exp
    return chi_square_sf_1df(chi2)


@dataclass
class GenotypeQC:
    maf_threshold: float = 0.01
    call_rate_threshold: float = 0.95
    hwe_p_threshold: float = 1e-6

    def compute_metrics(self, dataset: GenotypeDataset) -> Dict[str, Dict[str, float]]:
        metrics: Dict[str, Dict[str, float]] = {}
        maf = dataset.minor_allele_frequencies()
        call_rate = dataset.call_rate()
        heterozygosity = dataset.heterozygosity()
        for marker in dataset.markers:
            column = dataset._column(dataset.markers.index(marker))
            metrics[marker] = {
                "call_rate": call_rate[marker],
                "maf": maf[marker],
                "heterozygosity": heterozygosity[marker],
                "hwe_pvalue": hardy_weinberg_pvalue(column),
            }
        return metrics

    def filter(self, dataset: GenotypeDataset) -> Dict[str, object]:
        metrics = self.compute_metrics(dataset)
        fail_call = [marker for marker, values in metrics.items() if values["call_rate"] < self.call_rate_threshold]
        fail_maf = [marker for marker, values in metrics.items() if values["maf"] < self.maf_threshold]
        fail_hwe = [marker for marker, values in metrics.items() if values["hwe_pvalue"] < self.hwe_p_threshold]
        keep = [marker for marker in dataset.markers if marker not in set(fail_call + fail_maf + fail_hwe)]
        return {
            "metrics": metrics,
            "keep_markers": keep,
            "fail_call_rate": fail_call,
            "fail_maf": fail_maf,
            "fail_hwe": fail_hwe,
        }

    def apply(self, dataset: GenotypeDataset) -> GenotypeDataset:
        result = self.filter(dataset)
        keep_markers = result["keep_markers"]
        indices = [dataset.markers.index(marker) for marker in keep_markers]
        new_matrix = [[row[idx] for idx in indices] for row in dataset.matrix]
        return GenotypeDataset(list(dataset.individuals), keep_markers, new_matrix)
