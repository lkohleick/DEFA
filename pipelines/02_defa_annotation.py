#!/usr/bin/env python3
################################################
####  DEFA annotation: one study at a time  ####
################################################
"""Write a DEFA config per study and annotate each study with run.py, in sequence."""

import subprocess
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent))
from _shared import load_config, parse_args, run_script, selected_studies, study_paths

# the pipelines ship inside the repo, so run.py sits one level up
RUN_PY = Path(__file__).resolve().parent.parent / "run.py"


################################################
####               study config             ####
################################################

def dataset_block(cfg: dict, study: str, paths: dict) -> dict:
    """DEFA dataset paths of one study, from the config or from 01_ld_proxies output."""
    study_cfg = cfg["studies"][study]
    block     = {}

    # a study may bring its own lead and proxy tables instead of running 01
    if study_cfg.get("lead_variant_file") and study_cfg.get("proxy_file"):
        block["lead_variant_file"] = str(study_cfg["lead_variant_file"])
        block["proxy_file"]        = str(study_cfg["proxy_file"])
    elif paths["leads"].exists() and paths["proxies"].exists():
        block["lead_variant_file"] = str(paths["leads"])
        block["proxy_file"]        = str(paths["proxies"])
    else:
        raise ValueError(
            f"study {study} has no variants: set `lead_variant_file` plus `proxy_file`, "
            f"or run 01_ld_proxies.py first"
        )

    # the endpoint map of 01 is picked up when the study does not name its own
    endpoints = study_cfg.get("gwas_endpoint_file")
    if not endpoints and paths["variant_map"].exists():
        endpoints = paths["variant_map"]
    if endpoints:
        block["gwas_endpoint_file"] = str(endpoints)
    if study_cfg.get("true_positive_links_file"):
        block["true_positive_links_file"] = str(study_cfg["true_positive_links_file"])
    return block


def write_study_config(cfg: dict, study: str, paths: dict) -> Path:
    """Write the DEFA override config of one study and return its path."""
    study_cfg   = cfg["studies"][study]
    annotation  = cfg.get("annotation", {})
    run_graph   = annotation.get("run_graph_scorer", False)

    override = {
        "force_recompute": annotation.get("force_recompute", False),
        "task": {
            "name":             study,
            "output_dir":       str(paths["defa_dir"]),
            "run_graph_scorer": run_graph,
        },
        "dataset": dataset_block(cfg, study, paths),
    }

    graph_scorer = annotation.get("graph_scorer", {})
    if not isinstance(graph_scorer, dict):
        raise ValueError("annotation.graph_scorer is a block of DEFA settings, use run_graph_scorer to switch it on")
    if run_graph and graph_scorer:
        override["graph_scorer"] = graph_scorer

    # the two DEFA parameters a study may override, pipeline value first, study value wins
    parameters = {}
    threshold  = study_cfg.get("scaling_factor_threshold", annotation.get("scaling_factor_threshold"))
    if threshold is not None:
        parameters["scaling_factor_threshold"] = float(threshold)
    if study_cfg.get("interval_kb") is not None:
        parameters["interval_kb"] = int(study_cfg["interval_kb"])
    if parameters:
        override["parameters"] = parameters

    paths["defa_dir"].mkdir(parents=True, exist_ok=True)
    with open(paths["study_config"], "w") as fh:
        yaml.dump(override, fh, default_flow_style=False, sort_keys=False)
    return paths["study_config"]


################################################
####               annotation               ####
################################################

def run_defa(paths: dict, check_only: bool = False) -> None:
    """Run DEFA for one study, streaming its output to the console and the study log."""
    command = [
        sys.executable, str(RUN_PY),
        "--config", str(paths["study_config"]),
        "--cache-dir", str(paths["cache"]),
    ]
    # the preflight prints its own report, so only its verdict is passed on
    if check_only:
        command.append("--check-only")
        if subprocess.run(command).returncode != 0:
            raise SystemExit("\nfix the paths above, nothing was annotated")
        return

    # the run takes hours, so its output is written as it arrives
    with open(paths["defa_log"], "w") as log:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            print(line, end="")
            log.write(line)
    if process.wait() != 0:
        raise SystemExit(f"DEFA failed, see {paths['defa_log']}")


################################################
####                  main                  ####
################################################

def main() -> None:
    args = parse_args("Annotate each study's variants with DEFA")
    cfg  = load_config(args.config)

    studies = selected_studies(cfg, args.study)
    todo    = []

    # --> configs and preflight first, a bad path must not surface hours in <-- #
    for study in studies:
        paths = study_paths(cfg, study)
        if paths["defa_table"].exists() and not args.force:
            print(f"[{study}] variant table present, skipped (--force to redo)")
            continue
        write_study_config(cfg, study, paths)
        print(f"[{study}] checking configured paths", flush=True)
        run_defa(paths, check_only=True)
        todo.append((study, paths))

    # --> annotate, one study per process <-- #
    for study, paths in todo:
        print(f"\n[{study}] annotating", flush=True)
        run_defa(paths)
        print(f"[{study}] done: {paths['defa_table']}")


if __name__ == "__main__":
    run_script(main)
