import pandas as pd
from sqlmodel.ext.asyncio.session import AsyncSession
from sqlmodel import select
from sqlalchemy.orm import selectinload
from pathlib import Path
import uuid
import logging

from backend.models import models
from backend.analysis import qc, gp, gwas, mating, plots
from backend.schemas import schemas

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# In-memory cache for data and models
# For a production app, consider Redis or a more robust solution
CACHE = {
    "geno_matrix": None,
    "pheno_matrix": None,
    "marker_info": None,
    "grm": None,
    "gp_models": {}, # {model_id: model_object}
}

RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)


async def load_data_into_cache(db: AsyncSession):
    """Loads all relevant data from DB into pandas DataFrames in the cache."""
    logger.info("Loading data from database into cache...")

    # Genotypes
    stmt = select(models.Genotype)
    result = await db.execute(stmt)
    genotypes = result.scalars().all()
    if not genotypes:
        logger.warning("No genotype data found in the database.")
        return
    genotypes_df = pd.DataFrame([g.dict() for g in genotypes])

    # Phenotypes
    stmt = select(models.Phenotype)
    result = await db.execute(stmt)
    phenotypes = result.scalars().all()
    phenotypes_df = pd.DataFrame([p.dict() for p in phenotypes])

    # Markers
    stmt = select(models.Marker)
    result = await db.execute(stmt)
    markers = result.scalars().all()
    marker_info_df = pd.DataFrame([m.dict() for m in markers]).set_index('snp_id')

    # Pivot to wide format
    geno_matrix = genotypes_df.pivot(index='sample_id', columns='snp_id', values='gt')
    pheno_matrix = phenotypes_df.pivot(index='sample_id', columns='trait', values='value')

    CACHE["geno_matrix"] = geno_matrix
    CACHE["pheno_matrix"] = pheno_matrix
    CACHE["marker_info"] = marker_info_df
    logger.info(f"Data loaded. Genotypes: {geno_matrix.shape}, Phenotypes: {pheno_matrix.shape}")

async def run_qc_service(db: AsyncSession, params: schemas.QCParams):
    if CACHE["geno_matrix"] is None:
        await load_data_into_cache(db)
        if CACHE["geno_matrix"] is None:
            raise ValueError("Genotype data not available for QC.")

    geno_matrix_qc, summary, log_content = qc.run_qc_pipeline(CACHE["geno_matrix"], params.dict())

    # Update cache with QC'd data
    CACHE["geno_matrix_qc"] = geno_matrix_qc

    # Create a run record
    run = models.AnalyticsRun(run_type="QC", params=params.dict())
    db.add(run)
    await db.commit()
    await db.refresh(run)

    # Generate plots
    maf_series = qc.calculate_maf(geno_matrix_qc)
    maf_plot_fig = plots.plot_maf_distribution(maf_series)

    # Save plots and log
    run_dir = RESULTS_DIR / str(run.id)
    run_dir.mkdir(exist_ok=True)
    maf_plot_fig.savefig(run_dir / "maf_distribution.png")

    log_path = run_dir / "qc_log.txt"
    log_path.write_text(log_content)

    summary_with_log = {**summary, "log_file": str(log_path)}

    response = schemas.QCResponse(
        run_id=run.id,
        summary=schemas.QCResultSummary(**summary_with_log),
        plots={"maf_distribution": f"/api/results/{run.id}/plots/maf_distribution.png"}
    )
    return response

async def run_gp_service(db: AsyncSession, request: schemas.GPTrainRequest):
    if "geno_matrix_qc" not in CACHE:
        raise ValueError("Please run QC first.")

    geno_matrix = CACHE["geno_matrix_qc"]
    pheno_matrix = CACHE["pheno_matrix"]

    # For now, we handle one trait at a time as per the analysis function
    trait = request.traits[0]
    gp_request_dict = request.dict()
    gp_request_dict['traits'] = [trait]

    model, cv_results, gebvs = gp.run_gp_pipeline(geno_matrix, pheno_matrix, gp_request_dict)

    model_id = str(uuid.uuid4())
    CACHE["gp_models"][model_id] = {"model": model, "trait": trait}

    run = models.AnalyticsRun(run_type="GP", params=request.dict())
    db.add(run)
    await db.commit()
    await db.refresh(run)

    # Store CV results
    # In a real app, this would be a new table `ResultGpCv`

    # Store predicted GEBVs if any
    if gebvs:
        for gebv_data in gebvs:
            res = models.ResultGp(run_id=run.id, trait=trait, **gebv_data)
            db.add(res)
        await db.commit()

    return schemas.GPTrainResponse(
        run_id=run.id,
        model_id=model_id,
        cv_results={trait: cv_results}
    )

async def run_gwas_service(db: AsyncSession, request: schemas.GWASRunRequest):
    if "geno_matrix_qc" not in CACHE:
        raise ValueError("Please run QC first.")

    geno_matrix = CACHE["geno_matrix_qc"]
    pheno_matrix = CACHE["pheno_matrix"]
    marker_info = CACHE["marker_info"]

    gwas_results_df = gwas.run_gwas_pipeline(geno_matrix, pheno_matrix, request.dict())

    run = models.AnalyticsRun(run_type="GWAS", params=request.dict())
    db.add(run)
    await db.commit()
    await db.refresh(run)

    run_dir = RESULTS_DIR / str(run.id)
    run_dir.mkdir(exist_ok=True)

    # Save results table
    results_path = run_dir / "gwas_results.csv"
    gwas_results_df.to_csv(results_path)

    # Generate and save plots
    bonferroni_thresh = 0.05 / len(gwas_results_df)
    fdr_thresh = gwas_results_df[gwas_results_df['p_fdr'] < 0.05]['p_value'].max()

    manhattan_fig = plots.plot_manhattan(gwas_results_df, marker_info, bonferroni_thresh, fdr_thresh)
    qq_fig = plots.plot_qq(gwas_results_df)

    manhattan_fig.savefig(run_dir / "manhattan.png")
    qq_fig.savefig(run_dir / "qq.png")

    return schemas.GWASRunResponse(
        run_id=run.id,
        message="GWAS run completed.",
        result_table_path=f"/api/results/{run.id}/files/gwas_results.csv",
        plots={
            "manhattan_plot": f"/api/results/{run.id}/plots/manhattan.png",
            "qq_plot": f"/api/results/{run.id}/plots/qq.png"
        }
    )

async def run_mating_service(db: AsyncSession, request: schemas.MatingOptimizeRequest):
    if "geno_matrix_qc" not in CACHE:
        raise ValueError("Please run QC first.")

    geno_matrix = CACHE["geno_matrix_qc"]

    # Calculate GRM if not already in cache
    if CACHE.get("grm") is None:
        logger.info("Calculating GRM for mating optimization...")
        CACHE["grm"] = gp.calculate_grm(geno_matrix)
        logger.info("GRM calculation finished.")

    # We need GEBVs for all traits for the candidates
    # Let's assume they are already predicted and available in the pheno_matrix for simplicity
    # In a real scenario, you'd run GP for all traits on all candidates
    all_gebvs = CACHE["pheno_matrix"][list(request.trait_weights.keys())]

    top_pairs, plot_data = mating.run_mating_ga(CACHE["grm"], all_gebvs, request.dict())

    run = models.AnalyticsRun(run_type="MATING", params=request.dict())
    db.add(run)
    await db.commit()
    await db.refresh(run)

    # Save results
    for pair in top_pairs:
        res = models.ResultMating(run_id=run.id, **pair)
        db.add(res)
    await db.commit()

    return schemas.MatingOptimizeResponse(
        run_id=run.id,
        best_pairs=[schemas.MatingPair(**p) for p in top_pairs],
        plot_data=plot_data
    )
