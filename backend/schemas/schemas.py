from pydantic import BaseModel
from typing import List, Dict, Any, Optional

# --- Base Models ---
class SampleSchema(BaseModel):
    id: str

class PhenotypeSchema(BaseModel):
    sample_id: str
    trait: str
    value: float

class MarkerSchema(BaseModel):
    snp_id: str
    chr: str
    bp: int

class GenotypeSchema(BaseModel):
    sample_id: str
    snp_id: str
    gt: float

# --- API Request/Response Schemas ---

# /upload
class UploadResponse(BaseModel):
    message: str
    filename: str
    num_samples: int
    num_snps: int
    num_phenotypes: int

# /qc/run
class QCParams(BaseModel):
    max_sample_missing: float = 0.2
    max_snp_missing: float = 0.2
    min_maf: float = 0.05
    min_hwe_p: float = 1e-6

class QCResultSummary(BaseModel):
    initial_samples: int
    initial_snps: int
    samples_after_qc: int
    snps_after_qc: int
    mean_maf: float
    log_file: str

class QCResponse(BaseModel):
    run_id: int
    summary: QCResultSummary
    plots: Dict[str, str] # e.g., {"maf_distribution": "/api/results/plots/run_id/maf_dist.png"}

# /gp/train
class GPTrainRequest(BaseModel):
    traits: List[str]
    train_ids: List[str]
    candidate_ids: Optional[List[str]] = None
    k_folds: int = 5

class GPTrainResponse(BaseModel):
    run_id: int
    model_id: str
    cv_results: Dict[str, Dict[str, float]] # {trait: {r2, pearson_r, rmse}}

# /gp/predict
class GPPredictRequest(BaseModel):
    model_id: str
    predict_ids: List[str]

class GEBVSchema(BaseModel):
    sample_id: str
    gebv: float

class GPPredictResponse(BaseModel):
    run_id: int
    trait: str
    gebvs: List[GEBVSchema]

# /gwas/run
class GWASRunRequest(BaseModel):
    trait: str
    model: str = 'additive' # 'additive', 'dominant', 'recessive'
    method: str = 'ols' # 'ols', 'lmm'
    n_pcs: int = 5

class GWASRunResponse(BaseModel):
    run_id: int
    message: str
    result_table_path: str # url to download csv
    plots: Dict[str, str] # manhattan_plot, qq_plot

# /mating/optimize
class MatingOptimizeRequest(BaseModel):
    top_k: int = 20
    alpha: float = 1.0 # weight for performance
    beta: float = 1.0 # weight for diversity
    trait_weights: Dict[str, float]

class MatingPair(BaseModel):
    parent1: str
    parent2: str
    predicted_gain: float
    diversity_score: float

class MatingOptimizeResponse(BaseModel):
    run_id: int
    best_pairs: List[MatingPair]
    plot_data: Dict[str, List[float]] # { "gains": [...], "diversities": [...] }

# General Error
class ErrorDetail(BaseModel):
    code: str
    message: str
    details: Optional[Any] = None

class ErrorResponse(BaseModel):
    error: ErrorDetail
