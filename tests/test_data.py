import pytest

from breeder.data import GenotypeDataset, PhenotypeDataset


def test_genotype_loading_and_metrics():
    dataset = GenotypeDataset.from_csv("tests/data/genotype.csv", index_col="id")
    freqs = dataset.allele_frequencies()
    assert pytest.approx(freqs["M1"], rel=1e-9) == 0.4
    assert pytest.approx(freqs["M4"], rel=1e-9) == 0.4

    maf = dataset.minor_allele_frequencies()
    assert pytest.approx(maf["M2"], rel=1e-9) == 0.5

    call_rate = dataset.call_rate()
    assert pytest.approx(call_rate["M2"], rel=1e-9) == 0.8

    heterozygosity = dataset.heterozygosity()
    assert pytest.approx(heterozygosity["M1"], rel=1e-9) == 0.4


def test_phenotype_normalization_and_summary():
    phenotype = PhenotypeDataset.from_csv("tests/data/phenotype.csv", trait="yield", index_col="id")
    summary = phenotype.summary()
    assert pytest.approx(summary["mean"], rel=1e-3) == 11.92
    normalized = phenotype.normalize()
    assert pytest.approx(sum(normalized.values) / len(normalized.values), abs=1e-9) == 0.0
    mean_sq = sum(value**2 for value in normalized.values) / (len(normalized.values) - 1)
    assert pytest.approx(mean_sq, abs=1e-9) == 1.0
