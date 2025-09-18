import datetime
from typing import Optional, Any
from sqlmodel import Field, SQLModel, Relationship, JSON, Column


class Sample(SQLModel, table=True):
    id: str = Field(primary_key=True)

    phenotypes: list["Phenotype"] = Relationship(back_populates="sample")
    genotypes: list["Genotype"] = Relationship(back_populates="sample")

class Phenotype(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    sample_id: str = Field(foreign_key="sample.id")
    trait: str
    value: float

    sample: Sample = Relationship(back_populates="phenotypes")

class Marker(SQLModel, table=True):
    snp_id: str = Field(primary_key=True)
    chr: str
    bp: int

    genotypes: list["Genotype"] = Relationship(back_populates="marker")

class Genotype(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    sample_id: str = Field(foreign_key="sample.id", index=True)
    snp_id: str = Field(foreign_key="marker.snp_id", index=True)
    gt: float

    sample: Sample = Relationship(back_populates="genotypes")
    marker: Marker = Relationship(back_populates="genotypes")

class AnalyticsRun(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    run_type: str # 'QC', 'GP', 'GWAS', 'MATING'
    params: dict = Field(sa_column=Column(JSON))
    created_at: datetime.datetime = Field(default_factory=datetime.datetime.utcnow)
    results: list["ResultGwas"] = Relationship(back_populates="run")
    results_gp: list["ResultGp"] = Relationship(back_populates="run")
    results_mating: list["ResultMating"] = Relationship(back_populates="run")


# Example of a result table, you would create one for each analysis type
class ResultGwas(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    run_id: int = Field(foreign_key="analyticsrun.id")
    snp_id: str
    p_value: float
    effect_size: Optional[float] = None

    run: AnalyticsRun = Relationship(back_populates="results")

class ResultGp(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    run_id: int = Field(foreign_key="analyticsrun.id")
    sample_id: str
    trait: str
    gebv: float

    run: AnalyticsRun = Relationship(back_populates="results_gp")

class ResultMating(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    run_id: int = Field(foreign_key="analyticsrun.id")
    parent1: str
    parent2: str
    predicted_gain: float
    diversity_score: float

    run: AnalyticsRun = Relationship(back_populates="results_mating")
