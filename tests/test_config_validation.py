################################################
####    Unit tests: config path preflight    ####
################################################
"""Tests for the preflight that checks every configured path before a run starts."""

import pytest
from easydict import EasyDict as edict

from defa import config_validation as cv


def _set(config, dotted, value):
    # setattr, not setdefault: easydict only registers the attribute through __setattr__
    node = config
    parts = dotted.split(".")
    for part in parts[:-1]:
        if getattr(node, part, None) is None:
            setattr(node, part, edict())
        node = getattr(node, part)
    setattr(node, parts[-1], value)


@pytest.fixture(autouse=True)
def _managers_on_path(tmp_path, monkeypatch):
    # stub every env manager, so the preflight never depends on the test machine
    bin_dir = tmp_path / "manager_bin"
    bin_dir.mkdir()
    for name in cv.ENV_MANAGERS:
        stub = bin_dir / name
        stub.write_text("#!/bin/sh\n")
        stub.chmod(0o755)
    monkeypatch.setenv("PATH", str(bin_dir))


def _satisfied_config(tmp_path, graph_scorer=False):
    """A config whose every required path exists under tmp_path."""
    config = edict()
    _set(config, "task.run_graph_scorer", graph_scorer)

    # the variant input
    lead = tmp_path / "leads.tsv"
    lead.write_text("locus_id\n")
    proxies = tmp_path / "proxies.tsv"
    proxies.write_text("locus_id\tproxy_id\n")
    _set(config, "dataset.lead_variant_file", str(lead))
    _set(config, "dataset.proxy_file", str(proxies))
    _set(config, "dataset.true_positive_links_file", None)
    _set(config, "dataset.gwas_endpoint_file", None)

    # build a target for every entry in the requirement list, by kind
    requirements = cv._always_required() + cv._conditional(config)
    for i, (key, kind, _label) in enumerate(requirements):
        target = tmp_path / f"req_{i}"
        if kind == "dir":
            target.mkdir(exist_ok=True)
            _set(config, key, str(target))
        elif kind == "glob":
            target.mkdir(exist_ok=True)
            (target / "part_1.tsv.gz").write_text("x")
            _set(config, key, str(target / "*.tsv.gz"))
        elif kind == "exec":
            target.write_text("#!/bin/sh\n")
            target.chmod(0o755)
            _set(config, key, str(target))
        else:
            target.write_text("x")
            _set(config, key, str(target))

    # the CADD and plugin trees are derived from vep.cache_dir, not configured
    cache_dir = tmp_path / "vep_cache"
    (cache_dir / cv.CADD_SUBPATH).mkdir(parents=True, exist_ok=True)
    (cache_dir / cv.CADD_SUBPATH / "whole_genome_SNVs.tsv.gz").write_text("x")
    (cache_dir / cv.CADD_SUBPATH / "gnomad.genomes.r3.0.indel.tsv.gz").write_text("x")
    (cache_dir / cv.PLUGIN_SUBPATH).mkdir(parents=True, exist_ok=True)
    _set(config, "vep.cache_dir", str(cache_dir))
    _set(config, "vep.env_manager", "micromamba")

    return config


################################################
####            the raised error            ####
################################################


def test_validate_lists_every_problem_at_once(tmp_path):
    config = _satisfied_config(tmp_path)
    _set(config, "paths.omim.genes_catalogue", str(tmp_path / "gone_a"))
    _set(config, "paths.eqtl.eqtlgen", str(tmp_path / "gone_b"))
    _set(config, "vep.fasta_file", str(tmp_path / "gone_c"))

    with pytest.raises(cv.ConfigPathError) as excinfo:
        cv.validate_config_paths(config, tmp_path)

    message = str(excinfo.value)
    assert "3 configured path(s)" in message
    for key in ("paths.omim.genes_catalogue", "paths.eqtl.eqtlgen", "vep.fasta_file"):
        assert key in message
