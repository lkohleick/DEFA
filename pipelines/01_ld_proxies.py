#!/usr/bin/env python3
################################################
####   LD proxies: 1000 Genomes EUR, hg38   ####
################################################
"""Resolve lead loci by position in the 1000 Genomes EUR hg38 panel, compute LD proxies, and write DEFA inputs."""

import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).parent))
from _shared import load_config, parse_args, run_script, selected_studies, study_paths


################################################
####                settings                ####
################################################

# the 1000 Genomes EUR panel of LD_PANEL.md section 2
PANEL_BUILD  = "hg38"
LD_STAT      = "r2-unphased"
LEAD_COLUMNS =["locus_id", f"chrom_{PANEL_BUILD}", f"pos_{PANEL_BUILD}"]


@dataclass
class Settings:
    """Panel and LD settings, read from the `ld_proxies` block once."""
    panel_dir:    Path
    plink2:       str
    ld_window_kb: int
    ld_window_r2: float
    threads:      int


def settings_from_config(cfg: dict) -> Settings:
    """Read the `ld_proxies` block."""
    block = cfg.get("ld_proxies", {})
    panel = block.get("panel")

    # guard: the panel is the one thing with no sensible default
    if not panel or not panel.get("dir"):
        raise ValueError("ld_proxies.panel needs `dir`, see pipelines/config.template.yaml")

    return Settings(
        panel_dir    = Path(panel["dir"]),
        plink2       = block.get("plink2_binary", "plink2"),
        ld_window_kb = block.get("ld_window_kb", 3000),
        ld_window_r2 = block.get("ld_window_r2", 0.5),
        threads      = int(cfg.get("threads", 4)),
    )


def check_lead_columns(columns: list, leads_file: str) -> None:
    """Stop unless the leads file has the hg38 chrom/pos columns of the panel."""
    if columns != LEAD_COLUMNS:
        raise ValueError(
            f"{leads_file} must have columns {LEAD_COLUMNS}, got {columns}\n"
            f"lift the leads to {PANEL_BUILD} first"
        )


################################################
####             panel helpers              ####
################################################

def chrom_to_idx(chrom: str) -> str:
    """Map a plain chromosome code to the pipeline index, X -> 23."""
    # guard: one naming style only
    if chrom.startswith("chr"):
        raise ValueError(f"chromosome codes must not carry a chr prefix, got '{chrom}'")
    return "23" if chrom == "X" else chrom


def idx_to_ensembl(chrom_idx: str) -> str:
    """Map a pipeline index to the Ensembl chromosome name, 23 -> X."""
    return "X" if chrom_idx == "23" else chrom_idx


def pvar_path(st: Settings, chrom_idx: str) -> Path:
    """Path of the panel pvar for one chromosome index."""
    return st.panel_dir / f"plink2_chr{idx_to_ensembl(chrom_idx)}.pvar"


def read_pvar(st: Settings, chrom_idx: str, positions=None) -> pl.DataFrame:
    """Read CHROM/POS/ID/REF/ALT from a panel pvar, optionally only at `positions`."""
    path = pvar_path(st, chrom_idx)
    scan = pl.scan_csv(
        path, separator="\t", comment_prefix="##",
        schema_overrides={"#CHROM": pl.Utf8, "POS": pl.Int64, "ID": pl.Utf8, "REF": pl.Utf8, "ALT": pl.Utf8},
    ).select(["#CHROM", "POS", "ID", "REF", "ALT"])

    # scan the codes of the whole file and the requested rows together
    codes_q = scan.select(pl.col("#CHROM").unique())
    rows_q  = scan if positions is None else scan.filter(pl.col("POS").is_in(list(positions)))
    codes, rows = pl.collect_all([codes_q, rows_q])

    # guard: the file holds only its own chromosome, X may be coded 23 by --output-chr 26
    expected = {"X", "23"} if chrom_idx == "23" else {chrom_idx}
    found    = set(codes["#CHROM"].to_list())
    if not found <= expected or any(c.startswith("chr") for c in found):
        raise ValueError(f"{path} carries chromosome codes {sorted(found)}, expected {sorted(expected)}")
    return rows


################################################
####            lead resolution             ####
################################################

def resolve_leads(st: Settings, leads_file: str, paths: dict) -> None:
    """Match leads to panel variants by position and write per-chromosome SNP lists."""
    # --> read leads, which must be in hg38 <-- #
    leads = pl.read_csv(leads_file, separator="\t", schema_overrides={"locus_id": pl.Utf8})
    check_lead_columns(leads.columns, leads_file)
    _, lead_chrom, lead_pos = LEAD_COLUMNS
    leads = leads.with_columns(pl.col(lead_chrom).cast(pl.Utf8), pl.col(lead_pos).cast(pl.Int64))

    dups = leads.filter(pl.col("locus_id").is_duplicated())["locus_id"].unique().to_list()
    if dups:
        raise ValueError(f"duplicate locus_id in {leads_file}: {sorted(dups)[:20]}")

    # one record per lead, filled in below
    records = [
        {
            "locus_id": row["locus_id"], "chrom_idx": chrom_to_idx(row[lead_chrom]),
            "pos_panel": row[lead_pos], "panel_id": None,
            "status": None, "ref": None, "alt": None,
        }
        for row in leads.iter_rows(named=True)
    ]

    # --> match on (chrom, pos) per chromosome <-- #
    for chrom_idx in sorted({rec["chrom_idx"] for rec in records}, key=lambda c: (len(c), c)):
        chrom_recs = [rec for rec in records if rec["chrom_idx"] == chrom_idx]
        # no pvar for this chromosome, e.g. Y or MT
        if not pvar_path(st, chrom_idx).exists():
            for rec in chrom_recs:
                rec["status"] = "not_in_panel"
            continue

        hits = read_pvar(st, chrom_idx, positions={rec["pos_panel"] for rec in chrom_recs})
        by_pos = {pos: grp for (pos,), grp in hits.group_by(["POS"])}
        for rec in chrom_recs:
            grp = by_pos.get(rec["pos_panel"])
            if grp is None:
                rec["status"] = "not_in_panel"
            elif grp.height > 1:
                rec["status"] = "ambiguous_position"
            else:
                rec["panel_id"], rec["ref"], rec["alt"] = grp.row(0)[2:5]
                rec["status"] = "matched"

    # --> write resolution table <-- #
    resolved = pl.DataFrame(
        records,
        schema={
            "locus_id": pl.Utf8, "chrom_idx": pl.Utf8, "pos_panel": pl.Int64, "panel_id": pl.Utf8,
            "status": pl.Utf8, "ref": pl.Utf8, "alt": pl.Utf8,
        },
    ).with_columns(
        pl.col("chrom_idx").replace("23", "X").alias("chrom"),
    ).select(["locus_id", "panel_id", "chrom", "pos_panel", "status", "ref", "alt"])
    paths["ld_dir"].mkdir(parents=True, exist_ok=True)
    resolved.write_csv(paths["resolved"], separator="\t")
    print(resolved.group_by("status").len().sort("status"))

    # --> write snp lists for matched leads only <-- #
    paths["snp_lists"].mkdir(parents=True, exist_ok=True)
    # a rerun with fewer leads must not keep the lists of the previous one
    for stale in paths["snp_lists"].glob("snp.list.*.txt"):
        stale.unlink()
    matched = [rec for rec in records if rec["status"] == "matched"]
    for chrom_idx in {rec["chrom_idx"] for rec in matched}:
        ids = sorted({rec["panel_id"] for rec in matched if rec["chrom_idx"] == chrom_idx})
        (paths["snp_lists"] / f"snp.list.{chrom_idx}.txt").write_text("\n".join(ids) + "\n")


def chromosomes(paths: dict) -> list:
    """Chromosome indices with a SNP list for one study."""
    return [p.stem.replace("snp.list.", "") for p in sorted(paths["snp_lists"].glob("snp.list.*.txt"))]


################################################
####             ld computation             ####
################################################

def run_plink2(st: Settings, chrom_idx: str, paths: dict) -> None:
    """Compute LD proxies of one chromosome's leads against the whole panel."""
    out_prefix = paths["vcor_dir"] / f"ld.proxies.{chrom_idx}"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    command = [
        st.plink2,
        "--pfile", str(st.panel_dir / f"plink2_chr{idx_to_ensembl(chrom_idx)}"),
        f"--{LD_STAT}",
        "--ld-window-kb", str(st.ld_window_kb),
        "--ld-window-r2", str(st.ld_window_r2),
        "--threads", str(st.threads),
        "--ld-snp-list", str(paths["snp_lists"] / f"snp.list.{chrom_idx}.txt"),
        "--out", str(out_prefix),
    ]

    print(f"  chr{idx_to_ensembl(chrom_idx)}", flush=True)
    subprocess.run(command, check=True)


################################################
####              proxy tables              ####
################################################

def postprocess(st: Settings, phenotype_idx: str, paths: dict) -> None:
    """Key vcor output on locus_id and write the lead, proxy and endpoint tables for DEFA."""
    b = PANEL_BUILD
    resolved = pl.read_csv(
        paths["resolved"], separator="\t",
        schema_overrides={
            "locus_id": pl.Utf8, "panel_id": pl.Utf8, "chrom": pl.Utf8,
            "pos_panel": pl.Int64, "ref": pl.Utf8, "alt": pl.Utf8,
        },
    )
    matched = resolved.filter(pl.col("status") == "matched")

    # --> leads table <-- #
    matched.select(
        "locus_id",
        ("chr" + pl.col("chrom")).alias(f"chr_{b}"),
        pl.col("pos_panel").alias(f"pos_{b}"),
        pl.col("ref").alias(f"ref_{b}"),
        pl.col("alt").alias(f"alt_{b}"),
        pl.col("chrom").alias("chr_ensembl"),
        "panel_id",
    ).write_csv(paths["leads"], separator="\t")

    # --> proxies table, one vcor at a time <-- #
    proxy_columns = [
        "locus_id", "proxy_id", "r2", f"chr_{b}", f"pos_{b}", f"ref_{b}", f"alt_{b}",
        "panel_id",
    ]
    lead_keys = matched.select(["panel_id", "locus_id"])
    with open(paths["proxies"], "w") as fh:
        fh.write("\t".join(proxy_columns) + "\n")
        for chrom_idx in chromosomes(paths):
            vcor = pl.read_csv(
                paths["vcor_dir"] / f"ld.proxies.{chrom_idx}.vcor",
                separator="\t", schema_overrides={"#CHROM_A": pl.Utf8, "CHROM_B": pl.Utf8},
            )
            # recover locus_id from the lead's panel ID
            chunk = vcor.select(
                pl.col("ID_A").alias("panel_id"),
                pl.col("ID_B").alias("proxy_id"),
                pl.col("POS_B").cast(pl.Int64).alias("pos_panel"),
                pl.col("UNPHASED_R2").alias("r2"),
            ).join(lead_keys, on="panel_id", how="inner")

            # alleles from the panel pvar, keyed on position and ID
            alleles = read_pvar(st, chrom_idx, positions=set(chunk["pos_panel"].to_list())).select(
                pl.col("POS").alias("pos_panel"),
                pl.col("ID").alias("proxy_id"),
                pl.col("REF").alias("ref"),
                pl.col("ALT").alias("alt"),
            )
            chunk = chunk.join(alleles, on=["pos_panel", "proxy_id"], how="left")

            chunk.select(
                "locus_id",
                "proxy_id",
                "r2",
                pl.lit(f"chr{idx_to_ensembl(chrom_idx)}").alias(f"chr_{b}"),
                pl.col("pos_panel").alias(f"pos_{b}"),
                pl.col("ref").alias(f"ref_{b}"),
                pl.col("alt").alias(f"alt_{b}"),
                "panel_id",
            ).write_csv(fh, separator="\t", include_header=False)

    # --> endpoint map for every input locus <-- #
    resolved.select("locus_id").with_columns(
        pl.lit(phenotype_idx).alias("phenotype_idx")
    ).write_csv(paths["variant_map"], separator="\t")


################################################
####                  main                  ####
################################################

def main() -> None:
    args = parse_args("Resolve leads in the 1000 Genomes EUR hg38 panel and compute LD proxies")
    cfg  = load_config(args.config)
    st   = settings_from_config(cfg)

    for study in selected_studies(cfg, args.study):
        study_cfg = cfg["studies"][study]
        paths     = study_paths(cfg, study)

        # a study that brings its own variant file has nothing to resolve
        if not study_cfg.get("leads_file"):
            print(f"\n[{study}] no leads_file, skipped")
            continue
        if paths["proxies"].exists() and not args.force:
            print(f"\n[{study}] proxies present, skipped (--force to redo)")
            continue

        print(f"\n[{study}] resolving leads")
        resolve_leads(st, study_cfg["leads_file"], paths)

        print(f"[{study}] computing LD")
        for chrom_idx in chromosomes(paths):
            run_plink2(st, chrom_idx, paths)

        print(f"[{study}] writing lead, proxy and endpoint tables")
        postprocess(st, study_cfg.get("phenotype_idx", study), paths)
        print(f"[{study}] done: {paths['proxies']}")


if __name__ == "__main__":
    run_script(main)
