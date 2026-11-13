################################################
####     Unit tests: the variant table      ####
################################################
"""Tests for the join that assembles the output table."""

from types import SimpleNamespace

import polars as pl

from defa import create_variant_table as cvt

VEP_COLUMNS = [
    "VEP_max", "VEP_sum",
    "CADD_phred_max", "CADD_phred_mean", "CADD_phred_sum",
    "CADD_raw_max", "CADD_raw_mean", "CADD_raw_sum",
    "VEP_weighted_mean", "CADD_phred_weighted_mean", "CADD_raw_weighted_mean",
]

MARKER = "chr1:100_A_G"


################################################
####         the minimal input set          ####
################################################


def _genes_in_region(genes=("ENSG1", "ENSG2")):
    """Two candidate genes at one locus, the skeleton every join hangs off."""
    return pl.DataFrame({
        "locus_id":                       ["L1"] * len(genes),
        "marker_name_hg38":               [MARKER] * len(genes),
        "ensembl_gene_id":                list(genes),
        "hgnc_symbol":                    [g.replace("ENSG", "SYM") for g in genes],
        "distance_to_footprint":          [float(i * 1000) for i in range(len(genes))],
        "distance_to_tss":                [float(i * 1000) for i in range(len(genes))],
        "weighted_distance_to_footprint": [float(i * 1000) for i in range(len(genes))],
        "weighted_distance_to_tss":       [float(i * 1000) for i in range(len(genes))],
    })


def _vep(genes=("ENSG1", "ENSG2")):
    frame = pl.DataFrame({
        "locus_id":        ["L1"] * len(genes),
        "ensembl_gene_id": list(genes),
    })
    return frame.with_columns([pl.lit(1.0).alias(col) for col in VEP_COLUMNS])


def _build(**overrides):
    """Run create_variant_table over the minimal input set, with named parts swapped."""
    parts = {
        "variant_data": pl.DataFrame({"locus_id": ["L1"], "marker_name_hg38": [MARKER]}),
        "vep_annotations": _vep(),
        "genes_in_region": _genes_in_region(),
        "eqtl_catalogue_annotation": pl.DataFrame({
            "locus_id": ["L1"],
            "gene_id":  ["ENSG1"],
            "Adipose Subcutaneous": [0.9],
        }),
        "catlas_annotation": pl.DataFrame({
            "locus_id":  ["L1"],
            "gene_id":   ["ENSG1"],
            "ABC_score": [0.5],
        }),
        "ensembl_gene_info_df": pl.DataFrame({"ensembl_gene_id": ["ENSG1", "ENSG2"]}),
        "omim_genes_annotation": pl.DataFrame({
            "ensembl_gene_id": ["ENSG1"],
            "omim_gene":       [1],
        }),
        "ot_drug_targets_annotation": pl.DataFrame({
            "ensembl_gene_id":                   ["ENSG1"],
            "open_targets_advanced_drug_target": [1],
        }),
        "flames_annotations": {},
        "config": SimpleNamespace(
            dataset=SimpleNamespace(true_positive_links_file=None)
        ),
    }
    parts.update(overrides)
    return cvt.create_variant_table(**parts)


################################################
####              cardinality               ####
################################################


def test_a_duplicated_evidence_row_does_not_duplicate_the_output():
    # the same (locus, gene) twice, which a left join would fan out to two rows
    duplicated = pl.DataFrame({
        "locus_id": ["L1", "L1"],
        "gene_id":  ["ENSG1", "ENSG1"],
        "Adipose Subcutaneous": [0.9, 0.4],
    })
    table = _build(eqtl_catalogue_annotation=duplicated)
    assert table.height == 2


################################################
####             the zero-fills             ####
################################################


def test_true_positive_is_one_for_listed_links_and_zero_elsewhere(tmp_path):
    true_positives = tmp_path / "true_positives.tsv"
    true_positives.write_text("locus_id\tensembl_gene_id\nL1\tENSG1\n")
    config = SimpleNamespace(
        dataset=SimpleNamespace(true_positive_links_file=str(true_positives))
    )

    table = _build(config=config).sort("ensembl_gene_id")

    # the 0 arrives from the blanket fill_null, well after the join
    assert table["true_positive"].to_list() == [1, 0]
    # no true-positive row may be lost to the dedupe that follows the join
    assert table["true_positive"].sum() == 1
