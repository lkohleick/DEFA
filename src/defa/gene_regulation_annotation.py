################################################
####   Gene regulation: cATLAS ABC links    ####
################################################
"""cATLAS cCRE overlaps for each variant, joined to ABC scores and mapped to Ensembl genes."""

import os
from glob import glob
from typing import Dict, List, Tuple
import polars as pl
import scipy.io as sio
from pybedtools import BedTool


################################################
####             catlas inputs              ####
################################################


def load_catlas_data(config):
    """cATLAS cCRE coordinates and cell types, checked against the matrix shape."""
    # read the cCRE by cell type matrix
    mat = sio.mmread(config.paths.gene_regulation.catlas.matrix)
    # cCRE coordinates with a row index
    cCREs = pl.read_csv(
        config.paths.gene_regulation.catlas.cre_bed,
        separator="\t",
        has_header=False,
        new_columns=["chr", "start", "end"],
    ).with_columns(pl.arange(0, pl.count()).alias("cre_idx"))

    # cell types with a column index
    celltypes = pl.read_csv(
        config.paths.gene_regulation.catlas.celltypes,
        separator="\t",
        has_header=False,
        new_columns=["celltype"],
    ).with_columns(pl.arange(0, pl.count()).alias("ct_idx"))

    assert celltypes.shape[0] == mat.shape[1] and cCREs.shape[0] == mat.shape[0]
    return cCREs, celltypes


def process_abc_data(config):
    """Concatenate the per-cell-type ABC files, keeping the top-scoring link per cCRE."""
    abc_files = glob(config.paths.gene_regulation.catlas.abc_data)
    # cell type from the file name, highest ABC_score per cCRE
    return pl.concat(
        [
            pl.read_csv(f, separator="\t")
            .with_columns(
                pl.lit(os.path.basename(f).replace(".tsv.gz", "")).alias("celltype")
            )
            .group_by("cCRE")
            .agg(pl.all().sort_by("ABC_score").last())
            .drop("Promoter")
            for f in abc_files
        ]
    )


################################################
####            variant overlap             ####
################################################


def overlap(row: Dict, cCREs: pl.DataFrame) -> pl.DataFrame | None:
    """cCREs containing one variant, repeated alongside the variant row, or None."""
    # cCREs that strictly contain the position
    hit = (
        cCREs.filter(
            (pl.col("chr") == row["chr_hg38"])
            & (pl.col("start") < row["pos_hg38"])
            & (pl.col("end") > row["pos_hg38"])
        )
        .select(["start", "end", "cre_idx"])
        .rename({"start": "cre_start", "end": "cre_end"})
    )
    # repeat the variant row once per hit
    if hit.shape[0] > 0:
        row_df = pl.DataFrame(row)
        row_df = row_df.select(pl.all().repeat_by(pl.lit(hit.shape[0])))
        return pl.concat([row_df, hit], how="horizontal")
    return None


def explode_list_columns(
    df: pl.DataFrame, columns_to_explode: List[str]
) -> pl.DataFrame:
    """Explode the listed list columns, pass the rest through."""
    expressions = [
        pl.col(col).list.explode() if col in columns_to_explode else pl.col(col)
        for col in df.columns
    ]
    return df.select(expressions)


def process_variant_data(
    variant_data: pl.DataFrame, cCREs: pl.DataFrame, abc: pl.DataFrame
) -> pl.DataFrame:
    """Range-join variants onto cCREs in polars and attach the ABC links."""

    # variants strictly inside a cCRE
    filtered_data = variant_data.join_where(
        cCREs,
        pl.col("chr_hg38") == pl.col("chr"),
        pl.col("pos_hg38") < pl.col("end"),
        pl.col("start") < pl.col("pos_hg38"),
    )
    # build the cCRE name as chr:start-end
    data_with_cCRE_name = filtered_data.with_columns(
        pl.concat_str(
            [
                pl.col("chr_hg38"),
                pl.lit(":"),
                pl.col("start").cast(pl.Utf8),
                pl.lit("-"),
                pl.col("end").cast(pl.Utf8),
            ]
        ).alias("cCRE_name")
    )

    # join the ABC links
    abc_joined_data = data_with_cCRE_name.join(
        abc, left_on="cCRE_name", right_on="cCRE"
    )

    # drop duplicates
    return abc_joined_data.unique()

def process_variant_data_bed_tools(
    variant_data: pl.DataFrame, cCREs: pl.DataFrame, abc: pl.DataFrame
) -> pl.DataFrame:
    """Intersect variants with cCREs through bedtools and attach the ABC links."""
    # --> variants as BED <-- #
    # BED is 0-based, so start is pos - 1
    variant_data_bed = pl.DataFrame(
        {
            "chrom": variant_data.get_column("chr_hg38"),
            "start": variant_data.get_column("pos_hg38") - 1,
            "end": variant_data.get_column("pos_hg38"),
            "name": variant_data.get_column("marker_name_hg38"),
            "score": 0,
        }
    )
    variant_data_bed = BedTool.from_dataframe(variant_data_bed.to_pandas()).sort()
    variant_data_bed = variant_data_bed.sort()
    cCREs_bed = BedTool.from_dataframe(cCREs.to_pandas()).sort()
    cCREs_bed = cCREs_bed.sort()

    # --> intersect <-- #
    # overlapping cCREs for each variant
    filtered_data_bed = variant_data_bed.intersect(cCREs_bed, wa=True, wb=True)
    filtered_data = pl.DataFrame(filtered_data_bed.to_dataframe())
    # rename the bedtools columns back
    filtered_data = filtered_data.with_columns(
        [
            pl.col("chrom").alias("chr_hg38"),
            pl.col("end").alias("pos_hg38"),
            pl.col("end").alias("start_hg38"),
            pl.col("end").alias("end_hg38"),
            pl.col("name").alias("marker_name_hg38"),
            pl.col("thickStart").alias("start"),
            pl.col("thickEnd").alias("end"),
            pl.col("itemRgb").alias("cre_idx"),
        ]
    )

    # build the cCRE name as chr:start-end
    data_with_cCRE_name = filtered_data.with_columns(
        pl.concat_str(
            [
                pl.col("chr_hg38"),
                pl.lit(":"),
                pl.col("start").cast(pl.Utf8),
                pl.lit("-"),
                pl.col("end").cast(pl.Utf8),
            ]
        ).alias("cCRE_name")
    )

    # --> ABC links <-- #
    abc_joined_data = data_with_cCRE_name.join(
        abc, left_on="cCRE_name", right_on="cCRE"
    )

    # drop duplicates and rejoin the variant columns
    abc_joined_data = abc_joined_data.unique()
    result = variant_data.join(
        abc_joined_data, left_on="marker_name_hg38", right_on="marker_name_hg38"
    )
    result = result[
        [
            "locus_id",
            "variant_id",
            "chr_hg38",
            "pos_hg38",
            "start_hg38",
            "end_hg38",
            "ref_hg38",
            "chr_hg38_ensembl",
            "celltype",
            "Gene Name",
            "ABC_score",
            "Distance",
            "cCRE_name",
        ]
    ]

    return result


################################################
####        gene regulation evidence        ####
################################################


def collect_gene_regulation_data(
    variant_data: pl.DataFrame,
    biotype: List[str],
    config,
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """cATLAS ABC evidence per variant, mapped from transcript to Ensembl gene id."""
    # load cCREs and ABC links
    cCREs, _ = load_catlas_data(config)
    abc = process_abc_data(config)

    # overlap the variants with cCREs
    variant_data_catlas_annotation = process_variant_data_bed_tools(
        variant_data[

                [
                    "locus_id",
                    "variant_id",
                    "chr_hg38",
                    "pos_hg38",
                    "start_hg38",
                    "end_hg38",
                    "ref_hg38",
                    "alt_hg38",
                    "marker_name_hg38",
                    "chr_hg38_ensembl",
                ]

        ],
        cCREs,
        abc,
    )

    # --> transcript to gene <-- #
    def map_enst_to_ensg(
        catlas_evidence: pl.DataFrame, enst_to_ensg_map: pl.DataFrame,
    ):
        """Split `Gene Name` into symbol and transcript, then join the Ensembl gene id."""
        # `Gene Name` is symbol:transcript, version suffixes dropped
        catlas_evidence = catlas_evidence.with_columns(
            [
                pl.col("Gene Name")
                .str.split(":")
                .list.get(0)
                .str.replace(r"\.\d+$", "")
                .alias("hgnc_symbol"),
                pl.col("Gene Name")
                .str.split(":")
                .list.get(1)
                .str.replace(r"\.\d+$", "")
                .alias("transcript_id"),
            ]
        ).drop("Gene Name")
        catlas_evidence_to_ensg = catlas_evidence.join(
            enst_to_ensg_map,
            on="transcript_id",
            how="inner",
        )

        return catlas_evidence_to_ensg.unique()

    enst_to_ensg_map = pl.read_parquet(
        config.paths.gene_information.enst_to_ensg_map
    )

    variant_data_catlas_annotation = map_enst_to_ensg(
        variant_data_catlas_annotation, enst_to_ensg_map
    )

    return variant_data_catlas_annotation
