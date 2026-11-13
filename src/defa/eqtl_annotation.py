################################################
####   eQTL: catalogue evidence per gene    ####
################################################
"""eQTL Catalogue max-PIP evidence, scaled per variant and aggregated per (locus, gene)."""

from typing import List, Tuple

import polars as pl


################################################
####                 query                  ####
################################################


def query_eqtl_data(
    preprocessed_eqtl_data: pl.DataFrame,
    variant_data: pl.DataFrame,
    variant_col: str = "marker_name_hg38",
) -> pl.DataFrame:
    """Catalogue rows whose marker name matches a locus variant, one per (variant, gene)."""
    # query on marker name
    query = preprocessed_eqtl_data.filter(
        pl.col("marker_name_hg38").is_in(variant_data.select(variant_col).to_series())
    )
    query = query.unique(subset=[variant_col, "gene_id"])
    # TODO check how much we loose due to false marker_name -<
    return query


################################################
####              aggregation               ####
################################################


def process_eqtl_data(
    preprocessed_eqtl_data: pl.DataFrame,
    variant_data: pl.DataFrame,
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """Scale each tissue PIP by the variant `scaling_factor`, then sum and weighted-mean per (locus, gene)."""
    # --> join catalogue onto locus variants <-- #
    eqtl = query_eqtl_data(preprocessed_eqtl_data, variant_data)
    variant_data = variant_data.select(
        ["locus_id", "marker_name_hg38", "scaling_factor"]
    )
    variant_data_eqtl = variant_data.join(
        eqtl,
        on="marker_name_hg38",
        how="inner",
    ).unique(subset=["locus_id", "marker_name_hg38", "gene_id"])

    # --> scale <-- #
    # scale every tissue column, keys and the factor itself stay as they are
    cols_to_exclude = ["locus_id", "marker_name_hg38", "gene_id", "scaling_factor"]

    eqtl = variant_data_eqtl.with_columns(
        pl.all()
        .exclude(cols_to_exclude)
        .mul(pl.col("scaling_factor"))
        .name.suffix("_scaled")
    )
    # keep keys and scaled columns only
    eqtl = eqtl.select(pl.col(cols_to_exclude), pl.col("^.*_scaled$"))
    scaled_cols: List[str] = [c for c in eqtl.columns if c.endswith("_scaled")]

    # --> aggregate per (locus, gene) <-- #
    # sum of each scaled column plus the weight sum
    agg_exprs: List[pl.Expr] = [pl.col("scaling_factor").sum().alias("weight_sum")]
    for c in scaled_cols:
        agg_exprs.append(pl.col(c).sum().alias(f"{c}_sum"))

    eqtl_agg = eqtl.group_by(["locus_id", "gene_id"]).agg(agg_exprs)

    # weighted mean, sum over weight sum
    eqtl_agg = eqtl_agg.with_columns(
        [
            (pl.col(f"{c}_sum") / pl.col("weight_sum")).alias(
                f"{c.removesuffix('_scaled')}_weighted_mean"
            )
            for c in scaled_cols
        ]
    )

    return eqtl_agg#.drop([f"{c}_sum" for c in scaled_cols])
