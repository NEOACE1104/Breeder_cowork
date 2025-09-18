from pydantic_settings import BaseSettings
import os

class Settings(BaseSettings):
    DATABASE_URL: str = "sqlite+aiosqlite:///./app.db"
    CORS_ORIGINS: list[str] = ["http://localhost:5173", "http://localhost:3000"]
    SECRET_KEY: str = "a-very-secret-key-that-you-should-change"
    LOG_LEVEL: str = "INFO"
    MAX_UPLOAD_MB: int = 200

    class Config:
        env_file = ".env"
        env_file_encoding = 'utf-8'

settings = Settings()
