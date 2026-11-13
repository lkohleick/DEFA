################################################
####         Unit tests: VEP module         ####
################################################
"""Tests for the VEP output parsing and the per-locus weighting."""

from pathlib import Path

import pandas as pd
import pytest

from defa import vep_query_module as vq


################################################
####               vep output               ####
################################################


def _stub_run_vep_every_variant(input_file, output_file, *args, **kwargs):
    """Write one VEP-like output row per input line, like VEP does for each strand line."""
    # one row per input line
    header = "#Uploaded_variation\tGene\tExtra"
    extra  = "IMPACT=HIGH;CADD_PHRED=12;CADD_RAW=0.5;"
    rows   = []
    for line in Path(input_file).read_text().splitlines():
        chrom, pos, _end, ref_alt, _strand = line.split()
        rows.append(f"{chrom}_{pos}_{ref_alt}\tENSG0001\t{extra}")
    Path(output_file).write_text("\n".join([header, *rows]) + "\n")


def test_shared_proxy_counts_at_every_locus(monkeypatch, tmp_path):
    """A proxy shared by two loci counts once at each, with that locus's weight."""
    # one variant in the LD cloud of two loci, at a different r2 for each
    monkeypatch.setattr(vq, "run_vep_within_environment", _stub_run_vep_every_variant)
    variants_df = pd.DataFrame(
        {
            "chr_hg38_ensembl": ["1", "1"],
            "pos_hg38":         [10, 10],
            "ref_hg38":         ["A", "A"],
            "alt_hg38":         ["T", "T"],
            "locus_id":         ["l1", "l2"],
            "scaling_factor":   [0.8, 0.4],
        }
    )
    genes_df = pd.DataFrame({"ensembl_gene_id": ["ENSG0001"]})

    gene_summary, _ = vq.process_variants_with_cmd_vep(
        variants_df=variants_df,
        genes_df=genes_df,
        weight_col="scaling_factor",
        use_weights=True,
        build="GRCh38",
        VEP_path="/bin/vep",
        VEP_fa_cache="/tmp/fasta",
        VEP_dir_cache="/tmp/cache",
        outdir=tmp_path,
    )

    # each locus keeps the variant once, with its own weight
    summary = gene_summary.set_index("locus_id")
    assert sorted(summary.index) == ["l1", "l2"]
    assert summary.loc["l1", "weight_sum"] == pytest.approx(0.8)
    assert summary.loc["l2", "weight_sum"] == pytest.approx(0.4)
    assert summary.loc["l1", "VEP_sum"] == pytest.approx(0.8)
    assert summary.loc["l2", "CADD_phred_mean"] == pytest.approx(12)
