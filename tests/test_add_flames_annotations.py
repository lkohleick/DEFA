################################################
####     Unit tests: FLAMES annotations     ####
################################################
"""Tests for the FLAMES interval annotations."""

import os
from types import SimpleNamespace

import polars as pl
import pytest

from defa import add_flames_annotations as afa


################################################
####          interval annotations          ####
################################################


def _build_config(abc_dir="abc"):
    """Config namespace with the ABC, HACER and promoter paths."""
    return SimpleNamespace(
        paths=SimpleNamespace(
            gene_regulation=SimpleNamespace(
                abc_ep_data=os.path.join(abc_dir, "ABC_EP"),
                abc_crispr_data=os.path.join(abc_dir, "ABC_CRISPR"),
                hacer=os.path.join(abc_dir, "HACER"),
                promoter=os.path.join(abc_dir, "promoters.parquet"),
            )
        )
    )


def _write_abc_gene(directory, frame):
    """Write one ABC per-gene file, named after its `ENSG` column."""
    directory.mkdir(parents=True, exist_ok=True)
    frame.write_parquet(directory / f"{frame['ENSG'][0]}.parquet")


def test_get_abc_ep_annotations_scaling(tmp_path):
    """One variant inside one ABC-EP interval gives sum, max and mean of its score."""
    config = _build_config(str(tmp_path))
    # one variant, one candidate gene
    variant_data = pl.DataFrame(
        {
            "locus_id": ["l1"],
            "chr_hg19": ["1"],
            "pos_hg19": [5],
        }
    )
    genes_in_region = pl.DataFrame(
        {"locus_id": ["l1"], "ensembl_gene_id": ["ENSG1"]}
    )

    # one ABC-EP interval covering the variant
    _write_abc_gene(
        tmp_path / "ABC_EP",
        pl.DataFrame(
            {
                "chr": ["1"],
                "start": [1],
                "end": [10],
                "ENSG": ["ENSG1"],
                "ABC.Score": [0.4],
            }
        ),
    )

    result = afa.get_ABC_EP_annotations(variant_data, genes_in_region, config)

    # aggregate columns and the single score
    assert set(result.columns) == {
        "locus_id",
        "ensembl_gene_id",
        "ABC_EP_sum",
        "ABC_EP_max",
        "ABC_EP_mean",
    }
    assert result["ABC_EP_sum"][0] == pytest.approx(0.4)
