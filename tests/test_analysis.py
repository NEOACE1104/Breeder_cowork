import csv
import pytest

from breeder.analysis import (
    AdmixtureAnalyzer,
    AssociationAnalyzer,
    DendrogramAnalyzer,
    PopulationStructureAnalyzer,
)
from breeder.data import GenotypeDataset, PhenotypeDataset


def test_population_structure_analyzer_performs_pca():
    genotype = GenotypeDataset.from_csv("tests/data/genotype.csv", index_col="id")
    analyzer = PopulationStructureAnalyzer(n_components=2)
    analyzer.fit(genotype)
    scores = analyzer.transform(genotype)
    assert len(scores) == 5
    assert len(scores[0]) == 2
    means = [sum(component[i] for component in scores) / len(scores) for i in range(2)]
    assert pytest.approx(means[0], abs=1e-6) == 0.0
    variance = analyzer.explained_variance_ratio
    assert len(variance) == 2
    assert sum(variance) <= 1.0 + 1e-6


def test_association_analyzer_runs_linear_models():
    genotype = GenotypeDataset.from_csv("tests/data/genotype.csv", index_col="id")
    phenotype = PhenotypeDataset.from_csv("tests/data/phenotype.csv", trait="yield", index_col="id")
    with open("tests/data/covariates.csv", "r", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        cov = {row[0]: [float(row[1]), float(row[2])] for row in reader}
    analyzer = AssociationAnalyzer(min_call_rate=0.9)
    results = analyzer.run(genotype, phenotype, covariates=cov)
    markers = [entry["marker"] for entry in results]
    assert "M2" not in markers
    assert "M3" in markers
    m3 = next(entry for entry in results if entry["marker"] == "M3")
    assert 0.0 <= m3["pvalue"] <= 1.0
    assert not any(entry["effect"] != entry["effect"] for entry in results)  # no NaNs


def test_dendrogram_analyzer_produces_hierarchy():
    genotype = GenotypeDataset.from_csv("tests/data/genotype.csv", index_col="id")
    analyzer = DendrogramAnalyzer()
    analyzer.fit(genotype)
    root = analyzer.get_root()
    assert sorted(root.members()) == sorted(genotype.individuals)
    assert len(analyzer.merges_) == len(genotype.individuals) - 1


def test_admixture_analyzer_generates_membership():
    genotype = GenotypeDataset.from_csv("tests/data/genotype.csv", index_col="id")
    analyzer = AdmixtureAnalyzer(n_populations=2, max_iter=50, tol=1e-5, random_state=42)
    analyzer.fit(genotype)
    q = analyzer.q_matrix_
    assert len(q) == len(genotype.individuals)
    for row in q:
        assert pytest.approx(sum(row), abs=1e-6) == 1.0
        assert all(value >= 0 for value in row)
