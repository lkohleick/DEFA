################################################
####   Input files: id columns and schema   ####
################################################
"""Read a variant or true-positive-links input file, checking its id columns first."""

import polars as pl
from loguru import logger


################################################
####              the schema                ####
################################################

# id columns each input file must carry, after the legacy rename
REQUIRED_COLUMNS = {
    "lead":          ["locus_id"],
    "proxy":         ["locus_id", "proxy_id"],
    "true_positive": ["locus_id", "ensembl_gene_id"],
}

# pre-rename id columns, old name -> new name
LEGACY_COLUMNS = {
    "lead":          {"rsid": "locus_id"},
    "proxy":         {"rsid_lead": "locus_id", "rsid_proxy": "proxy_id"},
    "true_positive": {"rsid": "locus_id"},
}

ID_COLUMNS = ("locus_id", "proxy_id")


class InputSchemaError(ValueError):
    """Raised when an input file lacks a required id column."""


def read_header(path) -> list:
    """Column names of a tab-separated file, gzipped or not."""
    return pl.read_csv(path, separator="\t", n_rows=0).columns


def legacy_renames(columns, kind) -> dict:
    """Old id columns present in `columns` whose new name is not already there."""
    return {
        old: new for old, new in LEGACY_COLUMNS[kind].items()
        if old in columns and new not in columns
    }


def missing_columns(columns, kind) -> list:
    """Required id columns absent from `columns`, after the legacy rename."""
    renamed = {legacy_renames(columns, kind).get(c, c) for c in columns}
    return [c for c in REQUIRED_COLUMNS[kind] if c not in renamed]


################################################
####               the reader               ####
################################################


def read_input_table(path, kind) -> pl.DataFrame:
    """Read one input file with its id columns renamed to the current schema and cast to String."""
    columns = read_header(path)
    missing = missing_columns(columns, kind)
    if missing:
        raise InputSchemaError(
            f"{kind} file {path} lacks column(s) {missing}, found {columns}"
        )

    renames = legacy_renames(columns, kind)
    if renames:
        logger.warning(f"{kind} file {path} uses old column names, renaming {renames}")

    targets = set(ID_COLUMNS) | set(renames)
    overrides = {c: pl.String for c in columns if c in targets}

    df = pl.read_csv(
        path, separator="\t", infer_schema_length=10000, schema_overrides=overrides
    )
    return df.rename(renames)
