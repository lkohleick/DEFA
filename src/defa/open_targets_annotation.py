################################################
####     Open Targets: drug target flag     ####
################################################
"""Flag candidate genes that are Open Targets approved or late-stage drug targets."""

import polars as pl


def create_approved_drug_targets_annotation(
    genes_in_region: pl.DataFrame, open_targets_data_drug_targets_path: str
) -> pl.DataFrame:
    """`open_targets_advanced_drug_target` (0/1) per candidate gene."""
    # the target list holds approved and advanced clinical targets only
    open_targets_data_drug_targets = pl.read_csv(
        open_targets_data_drug_targets_path, separator="\t"
    )
    # flag candidate genes in the target list
    open_targets_data_drug_targets = genes_in_region.with_columns(
        [
            pl.col("ensembl_gene_id")
            .is_in(open_targets_data_drug_targets["id"])
            .cast(pl.Int8)
            .alias("open_targets_advanced_drug_target"),
        ]
    )

    return open_targets_data_drug_targets.select(
        [
            "ensembl_gene_id",
            "hgnc_symbol",
            "open_targets_advanced_drug_target",
        ]
    )
