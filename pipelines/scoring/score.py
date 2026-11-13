#!/usr/bin/env python3
################################################
####  Gene scoring: models and locus rank   ####
################################################
"""Score a DEFA variant table with the XGBoost models and rank the genes of each locus."""

import argparse
import os
import sys
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import xgboost as xgb

# `data` sits next to this file
sys.path.insert(0, str(Path(__file__).parent))
from data import VariantDataset, VariantDatasetConfig

# column present only when DEFA graph scoring was run
_GRAPH_FEATURE_SENTINEL = "graph_distance"

# distance column and the column that breaks its ties
_DIST_PAIRS = [
    ("distance_to_footprint",          "distance_to_tss"),
    ("distance_to_tss",                "distance_to_footprint"),
    ("weighted_distance_to_footprint", "weighted_distance_to_tss"),
    ("weighted_distance_to_tss",       "weighted_distance_to_footprint"),
]

# model families averaged into DEFA-mean
_FAMILIES = ("with_graph_features", "without_graph_features")


################################################
####                 models                 ####
################################################

def load_models(models_dir: str, use_graph: bool, use_no_graph: bool) -> dict:
    """Return {model_stem: fitted_estimator} for the selected model families."""
    selected = {}
    for fname in sorted(os.listdir(models_dir)):
        if not fname.endswith(".json"):
            continue

        # skip the families not selected
        is_no_graph = "without_graph_features" in fname
        if is_no_graph and not use_no_graph:
            continue
        if not is_no_graph and not use_graph:
            continue

        # xgboost's own JSON format, a pickle is bound to the xgboost version that wrote it
        model = xgb.XGBClassifier()
        model.load_model(Path(models_dir) / fname)
        selected[Path(fname).stem] = model
        print(f"  Loaded model: {fname}")

    # guard: at least one model
    if not selected:
        raise RuntimeError(
            "No models loaded. Check --models-dir and --use-graph-models / "
            "--use-no-graph-models flags."
        )
    return selected


################################################
####                 ranking                ####
################################################

def _rel_min_max_inverted(df: pd.DataFrame, col: str) -> pd.Series:
    """Per-locus min-max normalisation inverted, lowest value -> 1.0."""
    hi = df.groupby("locus_id")[col].transform("max")
    lo = df.groupby("locus_id")[col].transform("min")
    denom = hi - lo
    return np.where(denom > 0, (hi - df[col]) / denom, 1.0)


def _rel_min_max(df: pd.DataFrame, col: str) -> pd.Series:
    """Per-locus min-max normalisation, highest value -> 1.0."""
    hi = df.groupby("locus_id")[col].transform("max")
    lo = df.groupby("locus_id")[col].transform("min")
    denom = hi - lo
    return np.where(denom > 0, (df[col] - lo) / denom, 1.0)


# applies to distance measures
def _apply_tiebreaker(
    df: pd.DataFrame,
    rel_col: str,
    tiebreak_col: str,
) -> pd.Series:
    """Leave one 1.0 per locus, the tied gene with the lowest `tiebreak_col`."""
    rel = df[rel_col].copy()
    tied = rel == 1.0

    # minimum tiebreak distance among tied genes per locus
    tiebreak_min = (
        df[tiebreak_col].where(tied).groupby(df["locus_id"]).transform("min")
    )
    is_winner = tied & (df[tiebreak_col] == tiebreak_min)

    rel_tb = _rel_min_max_inverted(df, tiebreak_col)

    return np.where(
        is_winner,
        1.0,
        np.where(tied, rel_tb * 0.9999, rel),
    )


def rank_genes(df: pd.DataFrame, model_score_cols: List[str]) -> pd.DataFrame:
    """Add the per-locus `rel_*` scores of the distances and the models."""
    # --> distances, closest gene -> 1.0, one winner per locus <-- #
    for dist_col, tiebreak_col in _DIST_PAIRS:
        if dist_col not in df.columns:
            continue
        rel_col = f"rel_{dist_col}"
        df[rel_col] = _rel_min_max_inverted(df, dist_col)
        if tiebreak_col in df.columns:
            df[rel_col] = _apply_tiebreaker(df, rel_col, tiebreak_col)

    # --> models, highest probability -> 1.0, ties kept <-- #
    for col in model_score_cols:
        df[f"rel_{col}"] = _rel_min_max(df, col)

    return df


def add_defa_mean(df: pd.DataFrame, model_score_cols: List[str]) -> pd.DataFrame:
    """Add DEFA-mean, `defa_mean_*`, the mean of the `rel_` domain-model scores per family."""
    for family in _FAMILIES:
        rel_cols = [f"rel_{c}" for c in model_score_cols if c.endswith(f"_informed_{family}")]
        if not rel_cols:
            continue

        # average over the domain models of the family, 1.0 = top gene of every model
        df[f"defa_mean_{family}"] = df[rel_cols].mean(axis=1)

    return df


################################################
####                 scoring                ####
################################################

def score_table(
    input_path: str,
    models_dir: str,
    output_dir: str,
    study_name: str,
    use_graph_models: bool = True,
    use_no_graph_models: bool = True,
) -> dict:
    """Score one DEFA table and return the scored and best-genes paths."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --> load and preprocess <-- #
    print(f"\n[1/4] Loading {input_path}")
    cfg = VariantDatasetConfig(
        path=str(input_path),
        categorical_columns=["chr", "marker_name_hg38", "locus_id", "ensembl_gene_id", "hgnc_symbol"],
        columns_to_drop=["locus_id_right", "ensembl_gene_id_right"],
        number_of_random_columns=5,
    )
    dataset = VariantDataset(cfg).load().preprocess().encode_categorical()
    df = dataset.data

    # --> resolve the model families from the graph features <-- #
    has_graph_features = _GRAPH_FEATURE_SENTINEL in df.columns
    if not has_graph_features:
        print(
            f"  '{_GRAPH_FEATURE_SENTINEL}' column not found -- "
            "graph models disabled, no-graph models enabled automatically."
        )
    use_graph    = use_graph_models    and has_graph_features
    use_no_graph = use_no_graph_models or not has_graph_features

    # --> load models <-- #
    print(f"\n[2/4] Loading models from {models_dir}")
    print(f"  Graph models: {'yes' if use_graph else 'no'}  |  "
          f"No-graph models: {'yes' if use_no_graph else 'no'}")
    models = load_models(models_dir, use_graph, use_no_graph)

    # --> predict <-- #
    print(f"\n[3/4] Predicting with {len(models)} model(s)")
    model_score_cols = []
    for name, model in models.items():
        df[name] = model.predict_proba(df[model.feature_names_in_])[:, 1]
        model_score_cols.append(name)
        print(f"  {name}: done")

    # decode categoricals back to strings before writing
    dataset.decode_categorical()

    # drop the columns only the models read
    df = df.drop(columns=[c for c in df.columns if c.startswith("random_col_")])

    # --> rank genes per locus <-- #
    print("\n[4/4] Ranking genes per locus")
    df = rank_genes(df, model_score_cols)

    # DEFA-mean of the domain models, appended last
    df = add_defa_mean(df, model_score_cols)

    # --> write full scored table <-- #
    scored_path = out_dir / f"scored_table_{study_name}.tsv"
    df.to_csv(scored_path, sep="\t", index=False)
    print(f"\nScored table written: {scored_path}  ({len(df)} rows)")

    # --> write best-genes table, top gene per method <-- #
    # the distance and model scores of `rank_genes`, then DEFA-mean
    rel_cols  = [f"rel_{d}" for d, _ in _DIST_PAIRS if d in df.columns] + [f"rel_{m}" for m in model_score_cols]
    mean_cols = [f"defa_mean_{f}" for f in _FAMILIES if f"defa_mean_{f}" in df.columns]
    best_rows = []
    for col in rel_cols + mean_cols:
        # DEFA-mean has no rel_ column, its top genes hold the locus maximum, ties kept
        is_mean = col in mean_cols
        is_top  = df[col] == (df.groupby("locus_id")[col].transform("max") if is_mean else 1.0)
        top = df.loc[is_top, ["locus_id", "ensembl_gene_id", "hgnc_symbol"]].copy()
        top["rel_score"]      = 1.0
        top["scoring_method"] = col

        # carry the raw score, the mean itself for DEFA-mean
        raw_col = col if is_mean else col[len("rel_"):]
        if raw_col in df.columns:
            top["raw_score"] = df.loc[top.index, raw_col].values
        best_rows.append(top)

    best_genes = pd.concat(best_rows, ignore_index=True)
    best_path = out_dir / f"best_genes_{study_name}.tsv"
    best_genes.to_csv(best_path, sep="\t", index=False)
    print(f"Best-genes table written: {best_path}  ({len(best_genes)} rows)")
    print(f"  Methods: {sorted(best_genes['scoring_method'].unique())}")

    return {"scored": scored_path, "best_genes": best_path}


################################################
####                   cli                  ####
################################################

def main() -> None:
    parser = argparse.ArgumentParser(description="Score DEFA tables with XGBoost models")
    parser.add_argument("--input",               required=True, help="Path to defa_variant_table.tsv")
    parser.add_argument("--models-dir",          required=True, help="Directory containing *.json model files")
    parser.add_argument("--output-dir",          required=True, help="Directory to write output files")
    parser.add_argument("--study-name",          required=True, help="Study identifier used in output filenames")
    parser.add_argument("--use-graph-models",    action="store_true",
                        help="Include models trained with graph features")
    parser.add_argument("--use-no-graph-models", action="store_true",
                        help="Include models trained without graph features")
    args = parser.parse_args()

    score_table(
        input_path          = args.input,
        models_dir          = args.models_dir,
        output_dir          = args.output_dir,
        study_name          = args.study_name,
        use_graph_models    = args.use_graph_models,
        use_no_graph_models = args.use_no_graph_models,
    )


if __name__ == "__main__":
    main()
