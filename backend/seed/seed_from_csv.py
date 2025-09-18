import pandas as pd
from sqlmodel.ext.asyncio.session import AsyncSession
from sqlalchemy.dialects.sqlite import insert as sqlite_upsert
from backend.models import models
import logging

logger = logging.getLogger(__name__)

async def seed_data(db: AsyncSession, csv_file_path_or_buffer):
    """
    Seeds the database from a wide-format CSV file.
    The CSV is expected to have sample IDs, phenotype columns, and SNP columns.
    """
    logger.info("Starting data seeding from CSV...")
    df = pd.read_csv(csv_file_path_or_buffer)

    # Identify column types
    phenotype_cols = [c for c in df.columns if c.lower().startswith('trait_')]
    snp_cols = [c for c in df.columns if c.lower().startswith('snp_')]
    sample_col = 'sample_id' # Assuming this is the name of the sample ID column

    if sample_col not in df.columns:
        # If 'sample_id' is not a column, assume the first column is the sample ID.
        sample_col = df.columns[0]
        df.rename(columns={sample_col: 'sample_id'}, inplace=True)

    logger.info(f"Found {len(df)} samples, {len(phenotype_cols)} phenotypes, and {len(snp_cols)} SNPs.")

    # --- Prepare data for bulk insertion ---

    # Samples
    samples_data = df[[sample_col]].rename(columns={sample_col: 'id'})

    # Phenotypes
    phenotypes_data = df[[sample_col] + phenotype_cols].melt(
        id_vars=[sample_col],
        value_vars=phenotype_cols,
        var_name='trait',
        value_name='value'
    )
    phenotypes_data.rename(columns={sample_col: 'sample_id'}, inplace=True)
    # Remove 'trait_' prefix
    phenotypes_data['trait'] = phenotypes_data['trait'].str.replace('trait_', '', case=False)

    # Markers
    markers_data = []
    for snp_id in snp_cols:
        # Assuming format is SNP_CHR_BP, e.g., "snp_1_12345"
        try:
            _, chr_val, bp_val = snp_id.split('_')
            markers_data.append({'snp_id': snp_id, 'chr': chr_val, 'bp': int(bp_val)})
        except ValueError:
            # Fallback if format is different
            markers_data.append({'snp_id': snp_id, 'chr': 'unknown', 'bp': 0})

    # Genotypes
    genotypes_data = df[[sample_col] + snp_cols].melt(
        id_vars=[sample_col],
        value_vars=snp_cols,
        var_name='snp_id',
        value_name='gt'
    )
    genotypes_data.rename(columns={sample_col: 'sample_id'}, inplace=True)

    # --- Execute bulk inserts ---
    # Using SQLAlchemy Core for performance with SQLite
    # The `sqlite_upsert` handles cases where we re-upload data.

    # Clear existing data
    logger.info("Clearing existing genomic data...")
    for table in [models.Genotype, models.Phenotype, models.Sample, models.Marker]:
        await db.execute(table.__table__.delete())

    logger.info("Inserting new data...")
    if not samples_data.empty:
        stmt = sqlite_upsert(models.Sample).values(samples_data.to_dict(orient='records'))
        stmt = stmt.on_conflict_do_nothing(index_elements=['id'])
        await db.execute(stmt)

    if not markers_data:
        raise ValueError("No SNP markers found in the CSV file.")

    stmt = sqlite_upsert(models.Marker).values(markers_data)
    stmt = stmt.on_conflict_do_nothing(index_elements=['snp_id'])
    await db.execute(stmt)

    if not phenotypes_data.empty:
        await db.run_sync(
            lambda session: session.bulk_insert_mappings(models.Phenotype, phenotypes_data.to_dict(orient='records'))
        )

    if not genotypes_data.empty:
         await db.run_sync(
            lambda session: session.bulk_insert_mappings(models.Genotype, genotypes_data.to_dict(orient='records'))
        )

    await db.commit()
    logger.info("Data seeding completed successfully.")

    return len(samples_data), len(markers_data), len(phenotypes_data)
