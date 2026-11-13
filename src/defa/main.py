################################################
####  Workflow: annotation steps in order   ####
################################################
"""Run every annotation step in order and assemble the variant table."""

import sys
from pathlib import Path

import polars as pl
from loguru import logger

from .add_flames_annotations import get_flames_annotations
from .create_variant_table import create_variant_table
from .eqtl_annotation import process_eqtl_data
from .gene_annotation import get_genes_in_region
from .gene_regulation_annotation import collect_gene_regulation_data
from .graph_scoring import apply_graph_scores
from .omim_annotation import create_omim_gene_annotation
from .open_targets_annotation import create_approved_drug_targets_annotation
from .utils import save_or_load_step, setup_logging
from .variant_data_processing import collect_variant_data
from .vep_query_module import get_VEP_annotations


# put src/ on the path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root / "src"))


################################################
####                workflow                ####
################################################


def run_workflow(config, force_recompute: bool = False):
    """Run every annotation step, cached per step, and return the results by name."""
    setup_logging()
    logger.info("Starting main function")
    logger.info("Processing Variant Data")

    # --> variants and candidate genes <-- #
    # load leads and proxies
    variant_data = save_or_load_step(
        "variant_data",
        collect_variant_data,
        config,
        force_recompute=force_recompute,
    )

    # collect the genes in each locus region
    logger.info("Collecting genes in region.")
    genes_in_region = save_or_load_step(
        "genes_in_region",
        get_genes_in_region,
        variant_data,
        config,
        force_recompute=force_recompute,
    )

    logger.info("Getting VEP Annotations")

    # --> variant effect <-- #
    # query VEP
    vep_annotations,variant_level_vep_annotations = save_or_load_step(
        "vep_annotations",
        get_VEP_annotations,
        variant_data,
        genes_in_region,
        config,
        force_recompute=force_recompute,
    )


    # --> eQTL <-- #
    # load the eQTL Catalogue max-PIP table, shipped in the Zenodo deposit
    logger.info("Loading eQTL Catalogue max-PIP table")
    preprocessed_eqtl_data_df = pl.read_parquet(config.paths.eqtl.eqtl_catalogue_max_pip)

    # match variants against the eQTL table
    logger.info("Querying eQTL data")
    eqtl_catalgoue_annotations = save_or_load_step(
        "eqtl_annotations_variants",
        process_eqtl_data,
        preprocessed_eqtl_data_df,
        variant_data,
    )

    # --> gene regulation <-- #
    logger.info("Adding Gene Regulation Annotations")
    # load the transcript to gene map
    ensembl_gene_info_df =  pl.read_parquet(
        config.paths.gene_information.enst_to_ensg_map
    )
    # cATLAS ABC links
    logger.info("- Adding CATLAS Annotations")
    gene_regulation_evidence = save_or_load_step(
        "gene_regulation_evidence",
        collect_gene_regulation_data,
        variant_data,
        ensembl_gene_info_df,
        config,
    )
    # FLAMES annotations
    logger.info("- Adding Annotations from FLAMES")
    flames_annotations = save_or_load_step(
        "flames_annotations",
        get_flames_annotations,
        variant_data,
        genes_in_region,
        config,
        force_recompute=force_recompute,
    )

    # --> gene level <-- #
    # OMIM disease genes
    logger.info("Adding Gene Level Annotations")
    logger.info("- Adding OMIM Gene Annotations")
    omim_genes_annotation = save_or_load_step(
        "omim_gene_annotation",
        create_omim_gene_annotation,
        genes_in_region,
        config.paths.omim.genes_catalogue,
        force_recompute=force_recompute,
    )
    logger.info("- Adding Open Targets Gene Annotations")

    # Open Targets drug targets
    ot_drug_targets_annotation = save_or_load_step(
        "ot_drug_targets_annotation",
        create_approved_drug_targets_annotation,
        genes_in_region,
        config.paths.open_targets.drug_targets,
        force_recompute=force_recompute,
    )
    # collect the step outputs by name
    results = {
        "variant_data": variant_data,
        "variant_vep_annotations": vep_annotations,
        "genes_in_region": genes_in_region,
        "eqtl_catalgoue_annotations": eqtl_catalgoue_annotations,
        "gene_regulation_evidence": gene_regulation_evidence,
        "flames_annotations": flames_annotations,
        "omim_genes_annotation": omim_genes_annotation,
        "ot_drug_targets_annotation": ot_drug_targets_annotation,
    }

    # --> variant table <-- #
    # join all evidence onto the candidate genes
    logger.info("Integrating with Variant Table")
    variant_table = save_or_load_step(
        "variant_table",
        create_variant_table,
        variant_data= variant_data,
        vep_annotations=vep_annotations,
        genes_in_region=genes_in_region,
        eqtl_catalogue_annotation=eqtl_catalgoue_annotations,
        catlas_annotation=gene_regulation_evidence,
        ensembl_gene_info_df=ensembl_gene_info_df,
        omim_genes_annotation=omim_genes_annotation,
        ot_drug_targets_annotation=ot_drug_targets_annotation,
        flames_annotations=flames_annotations,
        config=config,
    )
    results["variant_table"] = variant_table
    # add graph scores when the toggle is on
    if config.task.run_graph_scorer:
        logger.info("Setting up Graph Scorer")
        variant_table, graph_scorer = apply_graph_scores(variant_table, config)
        results["graph_scorer"] = graph_scorer
        results["variant_table"] = variant_table

    logger.info("Workflow completed!")

    # drop empty results
    non_empty_results = {k: v for k, v in results.items() if v is not None}

    return non_empty_results


if __name__ == "__main__":
    # scalene is a dev-only dependency (pip install -e ".[dev]")
    from scalene import scalene_profiler

    scalene_profiler.start()
    results = run_workflow()
    scalene_profiler.stop()

    # print a sample of each result
    for key, value in results.items():
        print(f"{key} sample:")
        print(value.head())
