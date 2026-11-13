################################################
####   Config preflight: check every path   ####
################################################
"""Check every reference path a run needs before the run starts."""

import os
import shutil
from glob import glob
from pathlib import Path

from loguru import logger

from .input_schema import missing_columns, read_header
from .vep_query_module import ENV_MANAGERS

# shipped placeholder prefix, worth its own message: the config was never filled in
PLACEHOLDER_PREFIXES = ("/path/to", "path/to")

# derived from vep.cache_dir in run_vep_within_environment, not separately configurable
CADD_SUBPATH = "vep_plugins/CADD/krishna.gs.washington.edu/download/CADD/v1.6/GRCh38"
PLUGIN_SUBPATH = ".vep/Plugins38"


class ConfigPathError(Exception):
    """Raised when a path the run needs is missing, unreadable or still a placeholder."""


def _get(config, dotted):
    """Value at a dotted config key, or None if any level is missing."""
    node = config
    for part in dotted.split("."):
        node = getattr(node, part, None)
        if node is None:
            return None
    return node


################################################
####          the requirement list          ####
################################################


def _always_required():
    return [
        ("paths.gene_information.enst_to_ensg_map", "file", "transcript to gene map"),
        ("paths.omim.genes_catalogue", "file", "OMIM mim2gene"),
        ("paths.open_targets.drug_targets", "file", "Open Targets drug targets"),
        ("paths.eqtl.eqtl_catalogue_max_pip", "file", "eQTL Catalogue max-PIP table"),
        ("paths.eqtl.gtex_caviar_flames", "dir", "GTEx CAVIAR eQTL per-gene files"),
        ("paths.eqtl.eqtl_catalogue_rqtl", "dir", "eQTL Catalogue rQTL per-gene files"),
        ("paths.eqtl.eqtl_catalogue_eqtl", "dir", "eQTL Catalogue eQTL per-gene files"),
        ("paths.eqtl.eqtlgen", "dir", "eQTLGen per-gene files"),
        ("paths.gene_regulation.promoter", "file", "promoter regions"),
        ("paths.gene_regulation.abc_ep_data", "dir", "ABC enhancer-promoter per-gene files"),
        ("paths.gene_regulation.abc_crispr_data", "dir", "ABC CRISPR per-gene files"),
        ("paths.gene_regulation.fantom5", "dir", "FANTOM5 per-gene files"),
        ("paths.gene_regulation.genehancer", "dir", "GeneHancer per-gene files"),
        ("paths.gene_regulation.roadmap", "dir", "Roadmap Epigenomics per-gene files"),
        ("paths.gene_regulation.epimap", "dir", "EpiMap per-gene files"),
        ("paths.gene_regulation.cicero", "dir", "Cicero per-gene files"),
        ("paths.gene_regulation.hacer", "dir", "HACER directory"),
        ("paths.gene_regulation.hic_annotation.javierre_hic_data", "dir", "Javierre PCHi-C per-gene files"),
        ("paths.gene_regulation.hic_annotation.jung_hic_data", "dir", "Jung PCHi-C per-gene files"),
        ("paths.gene_regulation.catlas.matrix", "file", "CATlas matrix"),
        ("paths.gene_regulation.catlas.cre_bed", "file", "CATlas cCREs"),
        ("paths.gene_regulation.catlas.celltypes", "file", "CATlas cell types"),
        ("paths.gene_regulation.catlas.abc_data", "glob", "CATlas ABC scores"),
        ("vep.executable", "exec", "VEP executable"),
        ("vep.fasta_file", "file", "VEP reference FASTA"),
        ("vep.cache_dir", "dir", "VEP offline cache"),
        ("vep.env", "dir", "conda env holding VEP"),
    ]


def _conditional(config):
    """Requirements that depend on the task toggles."""
    reqs = []

    if bool(_get(config, "task.run_graph_scorer")):
        reqs += [
            ("graph_scorer.primekg_dir", "dir", "PrimeKG"),
            ("graph_scorer.txgnn_export_dir", "dir", "TxGNN export"),
            ("graph_scorer.ensembl_to_primekg_idx_file", "file", "gene to PrimeKG index map"),
            # graph_scoring reads this one unguarded
            ("dataset.gwas_endpoint_file", "file", "GWAS endpoint map"),
        ]
    elif _get(config, "dataset.gwas_endpoint_file"):
        reqs.append(("dataset.gwas_endpoint_file", "file", "GWAS endpoint map"))

    return reqs


################################################
####              the checks                ####
################################################


def _check_one(value, kind):
    """Return a problem string, or None when the path is usable."""
    if str(value).startswith(PLACEHOLDER_PREFIXES):
        return "still the shipped placeholder, fill it in"

    path = Path(value)
    if kind == "glob":
        return None if glob(str(value)) else "no files match this pattern"
    if not path.exists():
        return "does not exist"
    if kind == "dir" and not path.is_dir():
        return "is not a directory"
    if kind in ("file", "exec") and not path.is_file():
        return "is not a file"
    if kind == "exec" and not os.access(path, os.X_OK):
        return "is not executable"
    if kind in ("file", "exec") and not os.access(path, os.R_OK):
        return "is not readable"
    return None


def _check_schema(key, value, kind):
    """Return a problem string when the file lacks a required id column, else None."""
    try:
        columns = read_header(value)
    except Exception as e:
        return f"{key}  header could not be read: {e}\n    {value}"
    missing = missing_columns(columns, kind)
    if missing:
        return f"{key}  lacks column(s) {missing}, found {columns}\n    {value}"
    return None


def _check_input_files(inputs):
    """Existence, then id columns, of each (key, value, kind) input file."""
    problems = []
    for key, value, kind in inputs:
        problem = _check_one(value, "file")
        if problem:
            problems.append(f"{key}  {problem}\n    {value}")
            continue
        problem = _check_schema(key, value, kind)
        if problem:
            problems.append(problem)
    return problems


def _check_input_dataset(config):
    """The variant input: a lead file AND a proxy file."""
    lead = _get(config, "dataset.lead_variant_file")
    proxy = _get(config, "dataset.proxy_file")

    if proxy and not lead:
        return ["dataset.lead_variant_file  a proxy file was given without a lead file"]

    if lead and not proxy:
        return ["dataset.proxy_file  a lead file was given without a proxy file"]

    if not lead:
        return [
            "dataset  no variant input configured, set lead_variant_file plus proxy_file"
        ]

    return _check_input_files([
        ("dataset.lead_variant_file", lead,  "lead"),
        ("dataset.proxy_file",        proxy, "proxy"),
    ])


def _check_env_manager(config):
    """The tool that activates the VEP env: a known name, and on `PATH`."""
    manager = _get(config, "vep.env_manager")
    if manager not in ENV_MANAGERS:
        return [f"vep.env_manager  {manager!r} is not one of {sorted(ENV_MANAGERS)}"]
    if shutil.which(manager) is None:
        return [f"vep.env_manager  {manager} is not on PATH"]
    return []


def check_config_paths(config, project_root):
    """Every problem found, as a list of printable strings. Empty means good to run."""
    problems = _check_input_dataset(config)
    problems += _check_env_manager(config)

    # true-positive links are joined after every annotation, so check the header up front
    true_positives = _get(config, "dataset.true_positive_links_file")
    if true_positives:
        problems += _check_input_files(
            [("dataset.true_positive_links_file", true_positives, "true_positive")]
        )

    for key, kind, label in _always_required() + _conditional(config):
        value = _get(config, key)
        if value is None:
            problems.append(f"{key}  not set ({label})")
            continue
        problem = _check_one(value, kind)
        if problem:
            problems.append(f"{key}  {problem} ({label})\n    {value}")

    return problems


def validate_config_paths(config, project_root):
    """Raise ConfigPathError listing every unusable path, before the run starts."""
    problems = check_config_paths(config, project_root)

    # the CADD and plugin directories are derived from vep.cache_dir, not configured
    cache_dir = _get(config, "vep.cache_dir")
    if cache_dir and not str(cache_dir).startswith(PLACEHOLDER_PREFIXES) and Path(cache_dir).is_dir():
        cadd = Path(cache_dir) / CADD_SUBPATH
        for name in ("whole_genome_SNVs.tsv.gz", "gnomad.genomes.r3.0.indel.tsv.gz"):
            if not (cadd / name).exists():
                problems.append(f"vep CADD data  missing {name}\n    {cadd / name}")
        if not (Path(cache_dir) / PLUGIN_SUBPATH).is_dir():
            problems.append(f"vep plugins  directory missing\n    {Path(cache_dir) / PLUGIN_SUBPATH}")

    if problems:
        listing = "\n".join(f"  - {p}" for p in problems)
        raise ConfigPathError(
            f"{len(problems)} configured path(s) cannot be used:\n{listing}\n\n"
            "Fix these in your config, or pass --skip-path-check to run anyway."
        )

    logger.info("Config preflight passed, every configured path is present.")
