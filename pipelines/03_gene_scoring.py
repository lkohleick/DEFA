#!/usr/bin/env python3
################################################
####   Gene scoring: XGBoost model scoring  ####
################################################
"""Score each study's DEFA variant table and write the per-locus gene ranking."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent / "scoring"))
sys.path.append(str(Path(__file__).resolve().parent.parent / "src"))
from _shared import load_config, parse_args, run_script, selected_studies, study_paths
from scoring.score import score_table
from defa.utils import load_config as load_defa_config

# models of the DEFA data deposit, relative to `defa_data_dir`
DEFAULT_MODELS = "models/defa_models_v1.0.0"


################################################
####                  main                  ####
################################################

def main() -> None:
    args    = parse_args("Score DEFA variant tables with the trained XGBoost models")
    cfg     = load_config(args.config)
    scoring = cfg.get("scoring", {})

    # default to the deposit models under `defa_data_dir` of configs/base_config.local.yaml
    models_dir = scoring.get("models_dir") or str(Path(load_defa_config().defa_data_dir) / DEFAULT_MODELS)

    # guard: the models directory exists
    if not Path(models_dir).is_dir():
        raise ValueError(
            f"no models directory at {models_dir}, set `defa_data_dir` in "
            "configs/base_config.local.yaml or `scoring.models_dir` in the workflow config"
        )

    for study in selected_studies(cfg, args.study):
        paths = study_paths(cfg, study)

        if not paths["defa_table"].exists():
            raise SystemExit(f"[{study}] no variant table, run 02_defa_annotation.py first")
        if paths["scored"].exists() and not args.force:
            print(f"[{study}] scored table present, skipped (--force to redo)")
            continue

        print(f"\n[{study}] scoring {paths['defa_table']}")
        score_table(
            input_path          = paths["defa_table"],
            models_dir          = models_dir,
            output_dir          = str(paths["scoring_dir"]),
            study_name          = study,
            use_graph_models    = scoring.get("use_graph_models", True),
            use_no_graph_models = scoring.get("use_no_graph_models", True),
        )


if __name__ == "__main__":
    run_script(main)
