################################################
####   Scoring data: variant table loader   ####
################################################
"""Load a DEFA variant table and prepare it as model input for scoring."""

import functools
import time
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder


################################################
####           config and timing            ####
################################################


@dataclass
class VariantDatasetConfig:
    """Input path, categorical and dropped columns, and the number of extra model columns."""
    path: str
    categorical_columns: List[str] = field(default_factory=list)
    columns_to_drop: List[str] = field(default_factory=list)
    number_of_random_columns: int = 5


def timeit(func):
    """Print the wall time of each call."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        print(f"{func.__name__} took {time.time() - start:.2f}s")
        return result
    return wrapper


################################################
####                dataset                 ####
################################################


class VariantDataset:
    """Variant table with chained load, preprocess and categorical encoding steps."""
    def __init__(self, config: VariantDatasetConfig):
        self.config = config
        self._data: Optional[pd.DataFrame] = None
        self._encoders: dict = {}
        self._is_preprocessed = False
        self._is_encoded = False

    @timeit
    def load(self) -> "VariantDataset":
        """Read the tab-separated variant table."""
        if not __import__("pathlib").Path(self.config.path).exists():
            raise FileNotFoundError(f"Input file not found: {self.config.path}")
        self._data = pd.read_csv(self.config.path, sep="\t")
        return self

    @timeit
    def preprocess(self) -> "VariantDataset":
        """Add `chr`, drop merge leftovers and requested columns, add the model columns, fill nulls with 0."""
        if self._data is None:
            raise ValueError("Call load() first.")

        # chromosome number from marker_name_hg38, e.g. chr18:56911251_C_T, X as 23
        self._data["chr"] = (
            self._data["marker_name_hg38"]
            .str.split(":").str[0]
            .str.replace("chrX", "chr23")
            .str.extract(r"(\d+)")
            .astype(int)
        )

        # drop `_right` merge leftovers whose base column exists
        right_cols = [
            c for c in self._data.columns
            if c.endswith("_right") and c[:-6] in self._data.columns
        ]
        if right_cols:
            self._data.drop(columns=right_cols, inplace=True)

        # drop requested columns, skip missing ones
        to_drop = [c for c in self.config.columns_to_drop if c in self._data.columns]
        if to_drop:
            self._data.drop(columns=to_drop, inplace=True)

        rng = np.random.default_rng(seed=0)
        self._data["random_col_ones"] = 1
        for i in range(self.config.number_of_random_columns):
            self._data[f"random_col_{i}"] = rng.random(len(self._data))

        # fill nulls with 0
        self._data.fillna(0, inplace=True)
        self._is_preprocessed = True
        return self

    @timeit
    def encode_categorical(self, pre_existing_encoders: Optional[dict] = None) -> "VariantDataset":
        """Label-encode the categorical columns, reusing any encoder passed in."""
        if not self._is_preprocessed:
            raise ValueError("Call preprocess() first.")

        # cast the categorical columns to str
        self._encoders = dict(pre_existing_encoders or {})
        cols = [c for c in self.config.categorical_columns if c in self._data.columns]
        self._data[cols] = self._data[cols].astype(str)

        # transform with a known encoder, fit a new one otherwise
        for col in cols:
            if col in self._encoders:
                self._data[col] = self._encoders[col].transform(self._data[col])
            else:
                enc = LabelEncoder()
                self._data[col] = enc.fit_transform(self._data[col])
                self._encoders[col] = enc

        self._is_encoded = True
        return self

    @timeit
    def decode_categorical(self) -> "VariantDataset":
        """Map encoded columns back to their labels."""
        if not self._is_encoded:
            raise ValueError("Call encode_categorical() first.")
        for col, enc in self._encoders.items():
            if col in self._data.columns:
                self._data[col] = enc.inverse_transform(self._data[col])
        self._is_encoded = False
        return self

    @property
    def data(self) -> pd.DataFrame:
        """The current table."""
        return self._data
