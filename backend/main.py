from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from backend.core.config import settings
from backend.api.endpoints import router as api_router
from backend.db.database import init_db

app = FastAPI(
    title="Cucumber Genomic Breeding Web App",
    description="An application for genomic analysis of cucumber data.",
    version="0.1.0"
)

# Set up CORS
if settings.CORS_ORIGINS:
    app.add_middleware(
        CORSMiddleware,
        allow_origins=[str(origin) for origin in settings.CORS_ORIGINS],
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

@app.on_event("startup")
async def on_startup():
    """
    This function runs when the application starts.
    It initializes the database.
    """
    await init_db()

# Include the API router
app.include_router(api_router, prefix="/api")

@app.get("/")
async def root():
    return {"message": "Welcome to the Genomic Breeding API"}
