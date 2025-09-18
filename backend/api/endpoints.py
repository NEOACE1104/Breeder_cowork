from fastapi import APIRouter, Depends, UploadFile, File, HTTPException, status
from fastapi.responses import FileResponse
from sqlmodel.ext.asyncio.session import AsyncSession
import pandas as pd
from io import StringIO
import logging

from backend.db.database import get_async_session
from backend.schemas import schemas
from backend.services import analysis_service
from backend.seed.seed_from_csv import seed_data

router = APIRouter()
logger = logging.getLogger(__name__)

@router.post("/upload", response_model=schemas.UploadResponse)
async def upload_data(
    file: UploadFile = File(...),
    db: AsyncSession = Depends(get_async_session)
):
    """
    Uploads a CSV file, parses it, and seeds the database.
    The CSV should be an integrated format (sample, phenotype, snp1, snp2, ...).
    """
    if not file.filename.endswith('.csv'):
        raise HTTPException(status_code=400, detail="Invalid file type. Please upload a CSV.")

    try:
        contents = await file.read()
        csv_data = StringIO(contents.decode('utf-8'))

        num_samples, num_snps, num_phenotypes = await seed_data(db, csv_data)

        # After seeding, load data into cache for analysis
        await analysis_service.load_data_into_cache(db)

        return schemas.UploadResponse(
            message="File uploaded and database seeded successfully.",
            filename=file.filename,
            num_samples=num_samples,
            num_snps=num_snps,
            num_phenotypes=num_phenotypes
        )
    except Exception as e:
        logger.error(f"Error during file upload and seeding: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"An error occurred: {e}"
        )

@router.post("/qc/run", response_model=schemas.QCResponse)
async def run_qc(
    params: schemas.QCParams,
    db: AsyncSession = Depends(get_async_session)
):
    try:
        response = await analysis_service.run_qc_service(db, params)
        return response
    except Exception as e:
        logger.error(f"Error during QC run: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/gp/train", response_model=schemas.GPTrainResponse)
async def train_gp_model(
    request: schemas.GPTrainRequest,
    db: AsyncSession = Depends(get_async_session)
):
    try:
        response = await analysis_service.run_gp_service(db, request)
        return response
    except Exception as e:
        logger.error(f"Error during GP training: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/gwas/run", response_model=schemas.GWASRunResponse)
async def run_gwas_analysis(
    request: schemas.GWASRunRequest,
    db: AsyncSession = Depends(get_async_session)
):
    try:
        response = await analysis_service.run_gwas_service(db, request)
        return response
    except Exception as e:
        logger.error(f"Error during GWAS run: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/mating/optimize", response_model=schemas.MatingOptimizeResponse)
async def optimize_mating(
    request: schemas.MatingOptimizeRequest,
    db: AsyncSession = Depends(get_async_session)
):
    try:
        response = await analysis_service.run_mating_service(db, request)
        return response
    except Exception as e:
        logger.error(f"Error during Mating Optimization: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))

# Endpoint to serve result files (plots, csvs)
@router.get("/results/{run_id}/{file_type}/{filename}")
async def get_result_file(run_id: int, file_type: str, filename: str):
    """
    Serves a result file.
    file_type can be 'plots' or 'files'.
    """
    if file_type not in ["plots", "files"]:
        raise HTTPException(status_code=404, detail="File type not found.")

    file_path = analysis_service.RESULTS_DIR / str(run_id) / filename
    if not file_path.is_file():
        raise HTTPException(status_code=404, detail="File not found.")

    return FileResponse(file_path)
