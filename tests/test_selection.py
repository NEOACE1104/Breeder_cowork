import pytest

from breeder.data import GenotypeDataset, PhenotypeDataset
from breeder.selection import GenomicBLUP


def test_genomic_blup_fit_and_predict():
    genotype = GenotypeDataset.from_csv("tests/data/genotype.csv", index_col="id")
    phenotype = PhenotypeDataset.from_csv("tests/data/phenotype.csv", trait="yield", index_col="id")
    model = GenomicBLUP(alpha=0.1, random_state=0)
    model.fit(genotype, phenotype)
    preds = model.predict(genotype)
    assert len(preds) == len(genotype.individuals)
    mean_pred = sum(preds) / len(preds)
    assert mean_pred != pytest.approx(0.0, abs=1e-6)


def test_genomic_blup_cross_validation_metrics():
    genotype = GenotypeDataset.from_csv("tests/data/genotype.csv", index_col="id")
    phenotype = PhenotypeDataset.from_csv("tests/data/phenotype.csv", trait="yield", index_col="id")
    model = GenomicBLUP(alpha=0.1, random_state=0)
    cv = model.cross_validate(genotype, phenotype, n_splits=3, random_state=0)
    assert len(cv) == 3
    for entry in cv:
        assert set(entry.keys()) == {"fold", "rmse", "correlation", "n_test"}
        assert entry["rmse"] >= 0
