"""Command-line interface for the breeder toolkit."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Optional

from .analysis import (
    AdmixtureAnalyzer,
    AssociationAnalyzer,
    DendrogramAnalyzer,
    PopulationStructureAnalyzer,
)
from .data import GenotypeDataset, PhenotypeDataset
from .qc import GenotypeQC
from .selection import GenomicBLUP
from .plotting import admixture_plot, dendrogram_plot, gwas_manhattan, pca_scatter


def _load_genotype(path: str, index_col: Optional[str], sep: str) -> GenotypeDataset:
    return GenotypeDataset.from_csv(path, index_col=index_col, delimiter=sep)


def _load_phenotype(path: str, trait: str, index_col: Optional[str], sep: str) -> PhenotypeDataset:
    return PhenotypeDataset.from_csv(path, trait=trait, index_col=index_col, delimiter=sep)


def _load_covariates(path: str, index_col: Optional[str], sep: str) -> Dict[str, List[float]]:
    table: Dict[str, List[float]] = {}
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter=sep)
        header = next(reader)
        if index_col is None:
            index_idx = 0
        else:
            if index_col not in header:
                raise KeyError(f"Index column '{index_col}' not found in covariate file")
            index_idx = header.index(index_col)
        feature_indices = [i for i in range(len(header)) if i != index_idx]
        for row in reader:
            key = row[index_idx]
            values = [float(row[i]) for i in feature_indices]
            table[key] = values
    return table


def _write_csv(path: str, header: List[str], rows: List[List[object]]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def command_qc(args: argparse.Namespace) -> None:
    genotype = _load_genotype(args.genotype, args.geno_index, args.geno_sep)
    qc = GenotypeQC(
        maf_threshold=args.maf_threshold,
        call_rate_threshold=args.call_rate_threshold,
        hwe_p_threshold=args.hwe_p_threshold,
    )
    result = qc.filter(genotype)
    metrics_rows = [
        [marker, values["call_rate"], values["maf"], values["heterozygosity"], values["hwe_pvalue"]]
        for marker, values in result["metrics"].items()
    ]
    _write_csv(args.metrics_out, ["marker", "call_rate", "maf", "heterozygosity", "hwe_pvalue"], metrics_rows)
    if args.filtered_genotype:
        filtered = qc.apply(genotype)
        rows = [[ind] + [value if value is not None else "" for value in row] for ind, row in zip(filtered.individuals, filtered.matrix)]
        _write_csv(args.filtered_genotype, ["id"] + filtered.markers, rows)


def command_pca(args: argparse.Namespace) -> None:
    genotype = _load_genotype(args.genotype, args.geno_index, args.geno_sep)
    analyzer = PopulationStructureAnalyzer(n_components=args.components, scale=not args.no_scale)
    analyzer.fit(genotype)
    scores = analyzer.transform(genotype)
    header = ["individual"] + [f"PC{i+1}" for i in range(len(scores[0]))]
    rows = [[ind] + score for ind, score in zip(genotype.individuals, scores)]
    _write_csv(args.output, header, rows)
    if args.variance_out:
        variance_rows = [[f"PC{i+1}", ratio] for i, ratio in enumerate(analyzer.explained_variance_ratio)]
        _write_csv(args.variance_out, ["component", "explained_variance_ratio"], variance_rows)
    if args.plot_html:
        pca_scatter(scores, genotype.individuals, analyzer.explained_variance_ratio, args.plot_html)


def command_gwas(args: argparse.Namespace) -> None:
    genotype = _load_genotype(args.genotype, args.geno_index, args.geno_sep)
    phenotype = _load_phenotype(args.phenotype, args.trait, args.pheno_index, args.pheno_sep)
    covariates = _load_covariates(args.covariates, args.cov_index, args.cov_sep) if args.covariates else None
    analyzer = AssociationAnalyzer(min_call_rate=args.min_call_rate, impute=not args.no_impute)
    results = analyzer.run(genotype, phenotype, covariates=covariates)
    header = ["marker", "effect", "se", "t", "pvalue", "n"]
    rows = [[r["marker"], r["effect"], r["se"], r["t"], r["pvalue"], r["n"]] for r in results]
    _write_csv(args.output, header, rows)
    if args.plot_html:
        gwas_manhattan(results, args.plot_html)


def command_gblup(args: argparse.Namespace) -> None:
    genotype = _load_genotype(args.genotype, args.geno_index, args.geno_sep)
    phenotype = _load_phenotype(args.phenotype, args.trait, args.pheno_index, args.pheno_sep)
    model = GenomicBLUP(alpha=args.alpha, scale=not args.no_scale, n_splits=args.cv_splits, random_state=args.random_state)
    if args.cross_validate:
        if not args.output:
            raise SystemExit("--output is required for cross-validation results")
        results = model.cross_validate(genotype, phenotype, n_splits=args.cv_splits, random_state=args.random_state)
        header = ["fold", "rmse", "correlation", "n_test"]
        rows = [[r["fold"], r["rmse"], r["correlation"], r["n_test"]] for r in results]
        _write_csv(args.output, header, rows)
        return
    model.fit(genotype, phenotype)
    if args.output:
        predictions = model.predict(genotype)
        rows = [[ind, pred] for ind, pred in zip(genotype.individuals, predictions)]
        _write_csv(args.output, ["individual", "prediction"], rows)
    if args.predict_genotype:
        new_genotype = _load_genotype(args.predict_genotype, args.predict_index, args.geno_sep)
        predictions = model.predict(new_genotype)
        rows = [[ind, pred] for ind, pred in zip(new_genotype.individuals, predictions)]
        target = args.predict_output or "predictions.csv"
        _write_csv(target, ["individual", "prediction"], rows)


def command_dendrogram(args: argparse.Namespace) -> None:
    genotype = _load_genotype(args.genotype, args.geno_index, args.geno_sep)
    analyzer = DendrogramAnalyzer(linkage=args.linkage)
    analyzer.fit(genotype)
    rows = [
        [";".join(sorted(left)), ";".join(sorted(right)), distance]
        for left, right, distance in analyzer.merges_
    ]
    _write_csv(args.output, ["left", "right", "distance"], rows)
    if args.plot_html:
        dendrogram_plot(analyzer.get_root(), args.plot_html)


def command_admixture(args: argparse.Namespace) -> None:
    genotype = _load_genotype(args.genotype, args.geno_index, args.geno_sep)
    analyzer = AdmixtureAnalyzer(
        n_populations=args.populations,
        max_iter=args.max_iter,
        tol=args.tol,
        random_state=args.random_state,
    )
    analyzer.fit(genotype)
    header = ["individual"] + [f"cluster_{i+1}" for i in range(args.populations)]
    rows = [[ind] + list(props) for ind, props in zip(genotype.individuals, analyzer.q_matrix_)]
    _write_csv(args.output, header, rows)
    if args.allele_out:
        allele_header = ["marker"] + [f"cluster_{i+1}" for i in range(args.populations)]
        allele_rows = [
            [marker]
            + [analyzer.allele_frequencies_[k][idx] for k in range(args.populations)]
            for idx, marker in enumerate(genotype.markers)
        ]
        _write_csv(args.allele_out, allele_header, allele_rows)
    if args.plot_html:
        admixture_plot(analyzer.q_matrix_, genotype.individuals, args.plot_html)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="breeder", description="Scientific breeding analysis toolkit")
    subparsers = parser.add_subparsers(dest="command", required=True)

    qc_parser = subparsers.add_parser("qc", help="Compute genotype quality control metrics")
    qc_parser.add_argument("--genotype", required=True)
    qc_parser.add_argument("--geno-index", dest="geno_index", default=None)
    qc_parser.add_argument("--geno-sep", dest="geno_sep", default=",")
    qc_parser.add_argument("--metrics-out", dest="metrics_out", required=True)
    qc_parser.add_argument("--filtered-genotype", dest="filtered_genotype")
    qc_parser.add_argument("--maf-threshold", dest="maf_threshold", type=float, default=0.01)
    qc_parser.add_argument("--call-rate-threshold", dest="call_rate_threshold", type=float, default=0.95)
    qc_parser.add_argument("--hwe-p-threshold", dest="hwe_p_threshold", type=float, default=1e-6)
    qc_parser.set_defaults(func=command_qc)

    pca_parser = subparsers.add_parser("pca", help="Perform principal component analysis")
    pca_parser.add_argument("--genotype", required=True)
    pca_parser.add_argument("--geno-index", dest="geno_index", default=None)
    pca_parser.add_argument("--geno-sep", dest="geno_sep", default=",")
    pca_parser.add_argument("--components", type=int, default=10)
    pca_parser.add_argument("--output", required=True)
    pca_parser.add_argument("--variance-out", dest="variance_out")
    pca_parser.add_argument("--no-scale", action="store_true")
    pca_parser.add_argument("--plot-html", dest="plot_html")
    pca_parser.set_defaults(func=command_pca)

    gwas_parser = subparsers.add_parser("gwas", help="Run marker-wise association analysis")
    gwas_parser.add_argument("--genotype", required=True)
    gwas_parser.add_argument("--geno-index", dest="geno_index", default=None)
    gwas_parser.add_argument("--geno-sep", dest="geno_sep", default=",")
    gwas_parser.add_argument("--phenotype", required=True)
    gwas_parser.add_argument("--trait", required=True)
    gwas_parser.add_argument("--pheno-index", dest="pheno_index", default=None)
    gwas_parser.add_argument("--pheno-sep", dest="pheno_sep", default=",")
    gwas_parser.add_argument("--covariates")
    gwas_parser.add_argument("--cov-index", dest="cov_index", default=None)
    gwas_parser.add_argument("--cov-sep", dest="cov_sep", default=",")
    gwas_parser.add_argument("--min-call-rate", dest="min_call_rate", type=float, default=0.9)
    gwas_parser.add_argument("--no-impute", dest="no_impute", action="store_true")
    gwas_parser.add_argument("--output", required=True)
    gwas_parser.add_argument("--plot-html", dest="plot_html")
    gwas_parser.set_defaults(func=command_gwas)

    gblup_parser = subparsers.add_parser("gblup", help="Fit or evaluate a genomic BLUP model")
    gblup_parser.add_argument("--genotype", required=True)
    gblup_parser.add_argument("--geno-index", dest="geno_index", default=None)
    gblup_parser.add_argument("--geno-sep", dest="geno_sep", default=",")
    gblup_parser.add_argument("--phenotype", required=True)
    gblup_parser.add_argument("--trait", required=True)
    gblup_parser.add_argument("--pheno-index", dest="pheno_index", default=None)
    gblup_parser.add_argument("--pheno-sep", dest="pheno_sep", default=",")
    gblup_parser.add_argument("--alpha", type=float, default=1e-3)
    gblup_parser.add_argument("--no-scale", action="store_true")
    gblup_parser.add_argument("--cv-splits", type=int, default=5)
    gblup_parser.add_argument("--random-state", type=int, default=None)
    gblup_parser.add_argument("--cross-validate", action="store_true")
    gblup_parser.add_argument("--output")
    gblup_parser.add_argument("--predict-genotype", dest="predict_genotype")
    gblup_parser.add_argument("--predict-index", dest="predict_index", default=None)
    gblup_parser.add_argument("--predict-output", dest="predict_output", default="predictions.csv")
    gblup_parser.set_defaults(func=command_gblup)

    dendro_parser = subparsers.add_parser("dendrogram", help="Generate a dendrogram from genotype data")
    dendro_parser.add_argument("--genotype", required=True)
    dendro_parser.add_argument("--geno-index", dest="geno_index", default=None)
    dendro_parser.add_argument("--geno-sep", dest="geno_sep", default=",")
    dendro_parser.add_argument("--linkage", choices=["average", "single", "complete"], default="average")
    dendro_parser.add_argument("--output", required=True)
    dendro_parser.add_argument("--plot-html", dest="plot_html")
    dendro_parser.set_defaults(func=command_dendrogram)

    admixture_parser = subparsers.add_parser("admixture", help="Estimate admixture proportions")
    admixture_parser.add_argument("--genotype", required=True)
    admixture_parser.add_argument("--geno-index", dest="geno_index", default=None)
    admixture_parser.add_argument("--geno-sep", dest="geno_sep", default=",")
    admixture_parser.add_argument("--populations", type=int, required=True)
    admixture_parser.add_argument("--max-iter", dest="max_iter", type=int, default=500)
    admixture_parser.add_argument("--tol", type=float, default=1e-4)
    admixture_parser.add_argument("--random-state", dest="random_state", type=int, default=None)
    admixture_parser.add_argument("--output", required=True)
    admixture_parser.add_argument("--allele-out", dest="allele_out")
    admixture_parser.add_argument("--plot-html", dest="plot_html")
    admixture_parser.set_defaults(func=command_admixture)

    return parser


def main(argv: Optional[List[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
