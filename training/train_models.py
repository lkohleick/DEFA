################################################
####     Domain-informed model training     ####
################################################
"""Training procedure of the DEFA domain-informed XGBoost models"""

import random

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import average_precision_score
from sklearn.model_selection import RandomizedSearchCV
from sklearn.preprocessing import LabelEncoder

# label column, 1 = true effector gene of the locus
TARGET = "true_positive"
LOCUS  = "marker_name_hg38"
SEED   = 42

# label-encoded to integers before training
CATEGORICAL = ["chr", "marker_name_hg38", "locus_id", "ensembl_gene_id", "hgnc_symbol"]

# never used as features
EXCLUDE = ["locus_id", "ensembl_gene_id", "open_targets_advanced_drug_target"]

# additionally excluded for the `without_graph_features` models
GRAPH_FEATURES = [
    "graph_distance",
    "adamic_adar",
    "pagerank_score",
    "graph_distance_with_exclusions", # GO terms are excluded from the graph
    "adamic_adar_with_exclusions", # GO terms are excluded from the graph
    "pagerank_score_with_exclusions", # GO terms are excluded from the graph
    "txgnn_distmult",
]

# hyperparameter search, sampled per left-out chromosome
PARAM_GRID = {
    "learning_rate"   : [0.1],
    "max_depth"       : [3, 4, 5],
    "colsample_bytree": [0.8, 0.9],
    "min_child_weight": [3, 4, 5],
}
INT_PARAMS     = ["max_depth", "min_child_weight"]
N_ITER         = 10
CV             = 5
SCORING        = "accuracy"
EARLY_STOPPING = 10

# one locus in nine is held out from the final fit
HOLDOUT_SHARE = 9


################################################
####              preprocessing             ####
################################################

def prepare(table: pd.DataFrame) -> pd.DataFrame:
    """Turn a labelled DEFA variant table into an all-numeric training frame."""
    df = table.copy()

    # chromosome number from the marker name, chrX -> 23
    df["chr"] = df[LOCUS].str.split(":").str[0].str.replace("chrX", "chr23").str.extract(r"(\d+)")[0].astype(int)

    # keep only loci with at least one true positive gene
    positive_loci = df.loc[df[TARGET] == 1, LOCUS].unique()
    df = df[df[LOCUS].isin(positive_loci)]

    # drop merge leftovers whose base column exists
    df = df.drop(columns=[c for c in df.columns if c.endswith("_right") and c[:-6] in df.columns])

    rng = np.random.default_rng(SEED)
    df["random_col_ones"] = 1
    for i in range(5):
        df[f"random_col_{i}"] = rng.random(len(df))

    # missing evidence counts as zero
    df = df.fillna(0)

    # label-encode the id columns
    for col in CATEGORICAL:
        df[col] = LabelEncoder().fit_transform(df[col].astype(str))
    return df


################################################
####                training                ####
################################################

def tune_loco(data: pd.DataFrame, feats: list) -> list:
    """Tune one classifier per left-out chromosome and return each best parameter set."""
    best_params = []
    max_chrom   = data["chr"].max()
    for chrom in range(1, max_chrom + 1):
        # train on the loci of all other chromosomes
        keeps = set(data.loc[data["chr"] != chrom, LOCUS])

        # pad a small training set with loci off the mirrored chromosome
        if len(keeps) < 100:
            extra = set(data.loc[data["chr"] != (max_chrom + 1 - chrom), LOCUS])
            if len(extra) + len(keeps) >= 100:
                keeps.update(list(extra)[0 : 100 - len(keeps)])
            else:
                keeps.update(extra)

        # early-stop on the left-out loci
        train = data[data[LOCUS].isin(keeps)]
        valid = data[~data[LOCUS].isin(keeps)]

        # random search over PARAM_GRID
        clf = xgb.XGBClassifier(
            objective             = "binary:logistic",
            eval_metric           = "rmse",
            early_stopping_rounds = EARLY_STOPPING,
            seed                  = SEED,
            n_jobs                = -1,
        )
        search = RandomizedSearchCV(
            clf,
            param_distributions = PARAM_GRID,
            n_iter              = N_ITER,
            cv                  = CV,
            scoring             = SCORING,
            random_state        = SEED,
            n_jobs              = -1,
        )

        try:
            search.fit(train[feats], train[TARGET], eval_set=[(valid[feats], valid[TARGET])], verbose=False)
            best_params.append(search.best_params_)
        except Exception as e:
            print(f"Error training model for chromosome {chrom}: {e}")
    return best_params


def average_params(best_params: list) -> dict:
    """Average the per-chromosome parameter sets into one."""
    params = {}
    for key in PARAM_GRID:
        mean = np.mean([p[key] for p in best_params])
        # integers rounded, the rest to two decimals
        params[key] = int(np.round(mean)) if key in INT_PARAMS else round(float(mean), 2)
    return params


def fit_final(data: pd.DataFrame, feats: list, params: dict) -> tuple:
    """Fit the final model on eight of nine loci and score it on the ninth."""
    # shuffle loci and hold one in nine out
    loci = sorted(data[LOCUS].unique())
    random.seed(SEED)
    random.shuffle(loci)
    holdout = set(loci[: len(loci) // HOLDOUT_SHARE])
    train   = data[~data[LOCUS].isin(holdout)]
    valid   = data[data[LOCUS].isin(holdout)]

    # fit with the averaged parameters, early-stop on the holdout
    model = xgb.XGBClassifier(
        objective             = "binary:logistic",
        eval_metric           = "rmse",
        early_stopping_rounds = EARLY_STOPPING,
        seed                  = SEED,
        **params,
    )
    model.fit(
        train[feats], train[TARGET],
        eval_set = [(train[feats], train[TARGET]), (valid[feats], valid[TARGET])],
        verbose  = False,
    )

    # holdout average precision
    ap = average_precision_score(valid[TARGET], model.predict_proba(valid[feats])[:, 1])
    return model, ap


################################################
####                 entry                  ####
################################################

def train_domain_model(table: pd.DataFrame, with_graph_features: bool = True) -> dict:
    """Train one domain-informed model on a labelled DEFA variant table.

    """
    data = prepare(table)

    # every column is a feature except the excluded ones and the label
    exclude = EXCLUDE + ([] if with_graph_features else GRAPH_FEATURES) + [TARGET]
    feats   = [c for c in data.columns if c not in exclude]

    # tune per chromosome, average, refit once
    params    = average_params(tune_loco(data, feats))
    model, ap = fit_final(data, feats, params)
    print(f"params {params}, holdout average precision {ap:.3f}")
    return {"consensus_model": model}
