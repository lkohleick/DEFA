################################################
####     Command line: run the workflow     ####
################################################
"""Command-line entry point, parses the flags and writes the variant table."""

import argparse
import os
import sys
import shutil
from pathlib import Path

# put src/ on the path, run.py sits at the repo root
project_root = Path(__file__).parent
sys.path.append(str(project_root / "src"))

from defa.config_validation import ConfigPathError, validate_config_paths
from defa.main import run_workflow
from defa.utils import load_config, set_run_cache_dir


################################################
####                 run                    ####
################################################


def run():
    """Run the workflow on the global config and write the variant table."""
    global config
    # run all annotation steps
    results = run_workflow(config, force_recompute=config.force_recompute)

    name = config.task.name
    path = config.task.output_dir

    # create the output dir
    Path(path).mkdir(parents=True, exist_ok=True)

    # print a sample and write the variant table
    if "variant_table" in results:
        table = results["variant_table"]
        print("\nVariant table sample:")
        print(table.sample(min(5, table.height)))
        results["variant_table"].write_csv(
            f"{path}/defa_variant_table_{name}.tsv",
            separator="\t",
        )



################################################
####              entry point               ####
################################################


def main():
    """Parse the flags, run the preflight, clear caches, then run."""
    global config

    # --> flags <-- #
    parser = argparse.ArgumentParser(description="Run gene assignment workflow")
    parser.add_argument(
        "--config", help="Path to override configuration file (optional)"
    )
    parser.add_argument(
        "--force-recompute",
        action="store_true",
        help="Force recomputation of all steps",
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        help="Per-run cache directory for step pickles and VEP temp files. "
             "Use a distinct path for each parallel job (e.g. cache/run_job42). "
             "Defaults to cache/ when omitted.",
    )
    parser.add_argument(
        "--skip-path-check",
        action="store_true",
        help="Skip the preflight check that every configured reference path exists.",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Run the preflight path check and exit without annotating anything.",
    )

    args = parser.parse_args()
    # load the config with optional overrides
    config = load_config(args.config)

    # the flag overrides the config, it never turns recomputation off
    if args.force_recompute:
        config.force_recompute = True

    # resolve the per-run cache dir and register it with utils
    cache_dir = Path(args.cache_dir) if args.cache_dir else project_root / "cache"
    set_run_cache_dir(cache_dir)

    # --> preflight, before anything expensive starts <-- #
    if not args.skip_path_check:
        try:
            validate_config_paths(config, project_root)
        except ConfigPathError as e:
            print(f"\nConfig preflight failed.\n\n{e}\n", file=sys.stderr)
            sys.exit(2)
    if args.check_only:
        print("Preflight only, nothing was run.")
        sys.exit(0)

    # --> clear step caches <-- #
    # the shared graph scorer cache in the project-level cache/ is never deleted
    if config.force_recompute:
        CACHE_EXCLUSIONS = {
            "graph_scorer_graph.cache",
        }
        if cache_dir.exists():
            for item in os.listdir(cache_dir):
                if item in CACHE_EXCLUSIONS:
                    continue
                item_path = cache_dir / item
                if item_path.is_file():
                    item_path.unlink()
                elif item_path.is_dir():
                    shutil.rmtree(item_path)

    # run the workflow
    run()


if __name__ == "__main__":
    main()
