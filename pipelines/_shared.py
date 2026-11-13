################################################
####   Pipeline config, shared by 01 to 03  ####
################################################
"""Load the workflow config, reject unknown keys, and name every output path."""

import argparse
from pathlib import Path

import yaml


################################################
####               config keys              ####
################################################

# top level, one block per pipeline script
TOP_KEYS = {"output_dir", "threads", "studies", "ld_proxies", "annotation", "scoring"}

# per study, the input files and the two DEFA parameter overrides
STUDY_KEYS = {
    "leads_file", "phenotype_idx", "lead_variant_file", "proxy_file",
    "gwas_endpoint_file", "true_positive_links_file", "scaling_factor_threshold", "interval_kb",
}

# one block per script
LD_PROXIES_KEYS = {"panel", "plink2_binary", "ld_window_kb", "ld_window_r2"}
PANEL_KEYS      = {"dir"}
ANNOTATION_KEYS = {"run_graph_scorer", "graph_scorer", "force_recompute", "scaling_factor_threshold"}
SCORING_KEYS    = {"models_dir", "use_graph_models", "use_no_graph_models"}

TEMPLATE = "pipelines/config.template.yaml"


################################################
####              config loading            ####
################################################

def _unknown(block: dict, allowed: set, prefix: str) -> list:
    """Keys of one config block that no pipeline script reads."""
    return sorted(f"{prefix}{key}" for key in set(block or {}) - allowed)


def load_config(path: str) -> dict:
    """Read the workflow config and stop on a key no script reads."""
    with open(path) as fh:
        cfg = yaml.safe_load(fh) or {}

    # guard: required at the top level, every script needs both
    for key in ("output_dir", "studies"):
        if not cfg.get(key):
            raise ValueError(f"{path} must set `{key}`, see {TEMPLATE}")

    unknown = _unknown(cfg, TOP_KEYS, "")
    unknown += _unknown(cfg.get("ld_proxies"), LD_PROXIES_KEYS, "ld_proxies.")
    unknown += _unknown(cfg.get("ld_proxies", {}).get("panel"), PANEL_KEYS, "ld_proxies.panel.")
    unknown += _unknown(cfg.get("annotation"), ANNOTATION_KEYS, "annotation.")
    unknown += _unknown(cfg.get("scoring"), SCORING_KEYS, "scoring.")
    for study, study_cfg in cfg["studies"].items():
        unknown += _unknown(study_cfg, STUDY_KEYS, f"studies.{study}.")
    if unknown:
        raise ValueError(
            "config sets keys no pipeline script reads: " + ", ".join(unknown) + f"\nsee {TEMPLATE}"
        )

    # every path is resolved once, so a script may run from any directory
    cfg["output_dir"] = str(Path(cfg["output_dir"]).expanduser().resolve())
    return cfg


def study_paths(cfg: dict, study: str) -> dict:
    """Every input and output path of one study, the single source of the run layout."""
    root     = Path(cfg["output_dir"]) / study
    ld       = root / "ld_proxies"
    defa     = root / "defa"
    scoring  = root / "scoring"
    return {
        "ld_dir":       ld,
        "defa_dir":     defa,
        "scoring_dir":  scoring,
        # 01_ld_proxies
        "snp_lists":    ld / "tmp" / "snp_lists",
        "vcor_dir":     ld / "tmp" / "ld_proxies",
        "resolved":     ld / f"leads_resolved_{study}.tsv",
        "leads":        ld / f"leads_{study}.tsv",
        "proxies":      ld / f"proxies_{study}.tsv",
        "variant_map":  ld / f"variant_to_endpoint_map_{study}.tsv",
        # 02_defa_annotation
        "study_config": defa / "study_config.yaml",
        "cache":        defa / "cache",
        "defa_log":     defa / f"defa_{study}.log",
        "defa_table":   defa / f"defa_variant_table_{study}.tsv",
        # 03_gene_scoring
        "scored":       scoring / f"scored_table_{study}.tsv",
        "best_genes":   scoring / f"best_genes_{study}.tsv",
    }


################################################
####            command line                ####
################################################

def parse_args(description: str) -> argparse.Namespace:
    """Standard `--config`, `--study` and `--force` command line of all three scripts."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--config", required=True, help="Path to the workflow config YAML")
    parser.add_argument("--study", default=None, help="Run one study only (default: every study)")
    parser.add_argument(
        "--force", action="store_true", help="Redo studies whose outputs are already present"
    )
    return parser.parse_args()


def run_script(main) -> None:
    """Run one pipeline script, reporting a config mistake as a message rather than a traceback."""
    try:
        main()
    except ValueError as error:
        raise SystemExit(f"\nconfig error: {error}")


def selected_studies(cfg: dict, study: str | None) -> list:
    """Study names to run, in config order."""
    if study is None:
        return list(cfg["studies"])
    if study not in cfg["studies"]:
        raise ValueError(f"--study {study} is not in the config, have: {', '.join(cfg['studies'])}")
    return [study]
