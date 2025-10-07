from pathlib import Path

import pytest

from breeder.analysis import DendrogramAnalyzer
from breeder.data import GenotypeDataset
from breeder.plotting import admixture_plot, dendrogram_plot, gwas_manhattan, pca_scatter


@pytest.fixture
def genotype_dataset() -> GenotypeDataset:
    return GenotypeDataset.from_csv("tests/data/genotype.csv", index_col="id")


def test_pca_scatter_creates_html(tmp_path: Path, genotype_dataset: GenotypeDataset) -> None:
    from breeder.analysis import PopulationStructureAnalyzer

    analyzer = PopulationStructureAnalyzer(n_components=2)
    analyzer.fit(genotype_dataset)
    scores = analyzer.transform(genotype_dataset)
    target = tmp_path / "pca.html"
    pca_scatter(scores, genotype_dataset.individuals, analyzer.explained_variance_ratio, target)
    assert target.exists()
    content = target.read_text(encoding="utf-8")
    assert "Plotly.newPlot" in content


def test_gwas_manhattan_creates_html(tmp_path: Path) -> None:
    results = [
        {"marker": "M1", "pvalue": 0.05},
        {"marker": "M2", "pvalue": 1e-4},
        {"marker": "M3", "pvalue": 0.5},
    ]
    target = tmp_path / "gwas.html"
    gwas_manhattan(results, target)
    assert target.exists()
    assert "Plotly.newPlot" in target.read_text(encoding="utf-8")


def test_dendrogram_plot_creates_html(tmp_path: Path, genotype_dataset: GenotypeDataset) -> None:
    analyzer = DendrogramAnalyzer()
    analyzer.fit(genotype_dataset)
    target = tmp_path / "dendrogram.html"
    dendrogram_plot(analyzer.get_root(), target)
    assert target.exists()


def test_admixture_plot_creates_html(tmp_path: Path, genotype_dataset: GenotypeDataset) -> None:
    from breeder.analysis import AdmixtureAnalyzer

    analyzer = AdmixtureAnalyzer(n_populations=2, max_iter=10, random_state=0)
    analyzer.fit(genotype_dataset)
    target = tmp_path / "admixture.html"
    admixture_plot(analyzer.q_matrix_, genotype_dataset.individuals, target)
    assert target.exists()
