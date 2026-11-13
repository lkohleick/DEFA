################################################
####      Unit tests: workflow wiring       ####
################################################
"""Tests that `run_workflow` calls every step and returns the results by name."""

from types import SimpleNamespace

import polars as pl

from defa import main


################################################
####              config stub               ####
################################################


def _config(graph=False):
    """Minimal config namespace, graph scoring toggled by `graph`."""
    return SimpleNamespace(
        task=SimpleNamespace(run_graph_scorer=graph),
        paths=SimpleNamespace(
            eqtl=SimpleNamespace(
                eqtl_catalogue_max_pip="max_pip",
            ),
            gene_information=SimpleNamespace(enst_to_ensg_map="map"),
            omim=SimpleNamespace(genes_catalogue="omim"),
            open_targets=SimpleNamespace(drug_targets="ot"),
        ),
        parameters=SimpleNamespace(genome_build="hg38"),
        vep=SimpleNamespace(
            executable="/bin/vep",
            fasta_file="/tmp/fasta",
            cache_dir="/tmp/cache",
        ),
        graph_scorer=SimpleNamespace(
            primekg_dir="base",
            cache_file="cache",
            ensembl_to_primekg_idx_file="gene_idx",
            max_distance=4,
            alpha=0.85,
            exclude_node_types=None,
        ),
    )


################################################
####               step order               ####
################################################


def test_run_workflow_minimal_pipeline(monkeypatch):
    """Without graph scoring, every step result lands under its name."""
    config = _config()
    # cached output per step name
    stub_outputs = {
        "variant_data": "variants",
        "genes_in_region": "genes",
        "vep_annotations": ("gene_summary", "variant_level"),
        "eqtl_annotations_variants": "eqtl",
        "gene_regulation_evidence": "reg",
        "flames_annotations": "flames",
        "omim_gene_annotation": "omim",
        "ot_drug_targets_annotation": "ot",
        "variant_table": "table",
    }

    def fake_save_or_load_step(name, func, *args, **kwargs):
        """Return the stub for `name`, run `func` otherwise."""
        return stub_outputs.get(name, func(*args, **kwargs))

    # stub every step and the parquet reads
    monkeypatch.setattr(main, "save_or_load_step", fake_save_or_load_step)
    monkeypatch.setattr(main, "collect_variant_data", lambda *args, **kwargs: "variants")
    monkeypatch.setattr(main, "get_genes_in_region", lambda *args, **kwargs: "genes")
    monkeypatch.setattr(main, "get_VEP_annotations", lambda *args, **kwargs: ("gene_summary", "variant_level"))
    monkeypatch.setattr(main, "process_eqtl_data", lambda *args, **kwargs: "eqtl")
    monkeypatch.setattr(main, "collect_gene_regulation_data", lambda *args, **kwargs: "reg")
    monkeypatch.setattr(main, "get_flames_annotations", lambda *args, **kwargs: "flames")
    monkeypatch.setattr(main, "create_omim_gene_annotation", lambda *args, **kwargs: "omim")
    monkeypatch.setattr(main, "create_approved_drug_targets_annotation", lambda *args, **kwargs: "ot")
    monkeypatch.setattr(main, "create_variant_table", lambda *args, **kwargs: "table")
    monkeypatch.setattr(main.pl, "read_parquet", lambda *args, **kwargs: pl.DataFrame())

    # run and check the named results
    results = main.run_workflow(config, force_recompute=True)

    assert results["variant_data"] == "variants"
    assert results["variant_table"] == "table"
