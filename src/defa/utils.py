################################################
####   Utils: config, logging, step cache   ####
################################################
"""Config loading, logging setup and the per-step pickle cache."""

from pathlib import Path
from typing import Optional
import os
import pickle

import yaml
from easydict import EasyDict as edict
from loguru import logger

# per-run cache dir, set once by run.py
_run_cache_dir: Optional[Path] = None


################################################
####           logging and paths            ####
################################################


def set_run_cache_dir(path):
    """Register and create the per-run cache directory."""
    global _run_cache_dir
    _run_cache_dir = Path(path)
    _run_cache_dir.mkdir(parents=True, exist_ok=True)


def setup_logging():
    """Log to app.log and to stdout at INFO."""
    # drop the default handler
    logger.remove()
    # file handler
    logger.add("app.log", rotation="500 MB")
    # stdout handler
    logger.add(
        lambda msg: print(msg, end=""), colorize=True, level="INFO"
    )


def get_project_root():
    """Repo root, three levels above this file."""
    return Path(__file__).parent.parent.parent


################################################
####             config loading             ####
################################################


def deep_update(base_dict: dict, update_dict: dict) -> dict:
    """Recursively update `base_dict` in place with `update_dict`."""
    for key, value in update_dict.items():
        # recurse into nested dicts, overwrite everything else
        if (
            key in base_dict
            and isinstance(base_dict[key], dict)
            and isinstance(value, dict)
        ):
            deep_update(base_dict[key], value)
        else:
            base_dict[key] = value
    return base_dict


def unknown_keys(base_dict: dict, update_dict: dict, prefix: str = "") -> list:
    """Dotted paths an override sets that the base config does not define."""
    found = []
    for key, value in update_dict.items():
        dotted = f"{prefix}.{key}" if prefix else key
        # a key the base never declares is a typo or a name that has since been renamed
        if key not in base_dict:
            found.append(dotted)
        elif isinstance(base_dict[key], dict) and isinstance(value, dict):
            found.extend(unknown_keys(base_dict[key], value, dotted))
    return found


def merge_config_layer(config: dict, layer_path, base_config_path) -> dict:
    """Deep-merge one YAML layer, rejecting keys the base config does not define."""
    with open(layer_path, "r") as f:
        layer = yaml.safe_load(f)

    # reject keys the base does not know
    stale = unknown_keys(config, layer or {})
    if stale:
        listed = "\n".join(f"  {key}" for key in stale)
        raise ValueError(
            f"{layer_path} sets keys the base config does not define:\n"
            f"{listed}\n"
            f"check them against {base_config_path}"
        )

    return deep_update(config, layer or {})


# top-level keys whose value replaces `{key}` in every path of the merged config
PATH_ROOTS = ("defa_data_dir", "flames_data_dir")


def fill_path_roots(node: dict, roots: dict) -> dict:
    """Replace `{defa_data_dir}` and `{flames_data_dir}` in every string value."""
    for key, value in node.items():
        # recurse into nested dicts, substitute in strings
        if isinstance(value, dict):
            fill_path_roots(value, roots)
        elif isinstance(value, str):
            for name, root in roots.items():
                value = value.replace("{" + name + "}", root)
            node[key] = value
    return node


def load_config(override_config_path: Optional[str] = None) -> edict:
    """Load configs/base_config.yaml, deep-merged with the machine layer and an override."""
    base_config_path = get_project_root() / "configs" / "base_config.yaml"
    machine_config_path = get_project_root() / "configs" / "base_config.local.yaml"

    if not base_config_path.exists():
        raise FileNotFoundError(f"Base config not found at {base_config_path}")

    # read the base layer
    with open(base_config_path, "r") as f:
        config = yaml.safe_load(f)

    # the machine layer fills the shipped `/path/to` placeholders and is never committed
    if machine_config_path.exists():
        config = merge_config_layer(config, machine_config_path, base_config_path)

    # only the keys the override sets are replaced, the rest stay at their base values
    if override_config_path:
        config = merge_config_layer(config, override_config_path, base_config_path)

    # substitute the data roots into every path
    roots = {
        name: os.path.expanduser(str(config[name])).rstrip("/")
        for name in PATH_ROOTS
        if config.get(name)
    }
    config = fill_path_roots(config, roots)

    return edict(config)


################################################
####               step cache               ####
################################################


def get_cache_dir() -> Path:
    """Return the active run cache directory, creating it if needed."""
    if _run_cache_dir is not None:
        return _run_cache_dir
    # fall back to cache/ at the repo root
    cache_dir = get_project_root() / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def save_or_load_step(step_name, compute_func, *args, **kwargs):
    """Load `<step_name>.pkl` from the cache, or run `compute_func(*args, **kwargs)` and cache it."""
    filename = f"{step_name}.pkl"
    # pull the flag out so it never reaches compute_func
    force_recompute = kwargs.pop("force_recompute", False)

    result = save_or_load_pickle(
        lambda: compute_func(*args, **kwargs), filename, force_recompute=force_recompute
    )

    return result


def save_or_load_pickle(obj_or_func, filename, force_recompute=False):
    """Load `filename` from the cache, or pickle `obj_or_func` (called if callable) and return it."""
    pickle_file_path = get_cache_dir() / filename

    # load the cached pickle, recompute if it is truncated or corrupt
    if pickle_file_path.exists() and not force_recompute:
        try:
            with open(pickle_file_path, "rb") as f:
                return pickle.load(f)
        except (EOFError, pickle.UnpicklingError) as e:
            logger.warning(f"Error loading pickle file {filename}: {e}")
            logger.info("Recomputing and saving new pickle file.")
            force_recompute = True

    # compute or take the object as given
    if callable(obj_or_func):
        obj = obj_or_func()
    else:
        obj = obj_or_func

    try:
        with open(pickle_file_path, "wb") as f:
            pickle.dump(obj, f)
    except Exception as e:
        logger.error(f"Error saving pickle file {filename}: {e}")

    return obj
