################################################
####   Unit tests: pure transformations     ####
################################################
"""Tests for the functions that take a frame and return a frame."""

import polars as pl

from defa import gene_annotation as ga
from defa import variant_data_processing as vdp


################################################
####          chromosome handling           ####
################################################


def test_create_marker_name_orders_alleles_and_prefix():
    df = pl.DataFrame({"CHR": ["1"], "POS": [123], "REF": ["T"], "ALT": ["A"]})
    marker = vdp.create_marker_name(df, genome_build="hg38")
    # the alleles are sorted, so T/A comes back as A_T
    assert marker.to_list()[0] == "chr1:123_A_T"


################################################
####            region centroid             ####
################################################


def test_compute_regions_weights_the_centroid_by_scaling_factor():
    # positions 100 and 200 at weights 0.25 and 0.75 put the centre of gravity at 175,
    # not at the unweighted midpoint of 150
    df = pl.DataFrame(
        {
            "locus_id": ["L1", "L1"],
            "start_hg38": [100, 200],
            "end_hg38": [110, 210],
            "scaling_factor": [0.25, 0.75],
        }
    )
    out = ga.compute_regions(df)

    assert out.height == 1
    assert out["region_start"][0] == 100
    assert out["region_end"][0] == 210
    assert out["num_variants"][0] == 2
    assert out["region_centroid_weighted"][0] == 175
