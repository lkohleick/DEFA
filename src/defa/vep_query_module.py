################################################
####      VEP: variant effect per gene      ####
################################################
"""Run command-line VEP with CADD on both allele orientations and summarise per (locus, gene)."""

import os
import subprocess
import tempfile

import pandas as pd
import polars as pl
from loguru import logger

from .utils import get_cache_dir

# shell prefix that activates the VEP env, per package manager
# mamba 1.x has no `shell` command, so it activates through the conda it ships with
ENV_MANAGERS = {
    "micromamba": 'eval "$(micromamba shell hook -s bash)" && micromamba activate {env}',
    "mamba":      (
        'if mamba shell --help >/dev/null 2>&1; '
        'then eval "$(mamba shell hook -s bash)" && mamba activate {env}; '
        'else eval "$(conda shell.bash hook)" && conda activate {env}; fi'
    ),
    "conda":      'eval "$(conda shell.bash hook)" && conda activate {env}',
}


################################################
####            command-line vep            ####
################################################


def _map_to_loci(vep_df, variant_map, keep_weight):
    """Attach locus, variant id and weight to VEP rows, one copy per locus of the variant."""
    # join the map on the uploaded variation id
    cols   = ["vcf_id", "variant_id", "locus_id"] + (["scaling_factor"] if keep_weight else [])
    vep_df = vep_df.merge(
        variant_map[cols], left_on="#Uploaded_variation", right_on="vcf_id", how="left"
    ).drop(columns="vcf_id")

    # an uploaded variant missing from the map keeps the old default weight
    if keep_weight:
        vep_df["scaling_factor"] = vep_df["scaling_factor"].fillna(1.0)
    return vep_df


def process_variants_with_cmd_vep(
    variants_df,
    genes_df,
    weight_col=None,
    use_weights=False,
    build="GRCh38",
    VEP_path=None,
    VEP_fa_cache=None,
    VEP_dir_cache=None,
    VEP_env=None,
    VEP_env_manager="micromamba",
    VEP_fork=8,
    outdir="vep_output",
):
    """Return (gene summary, variant-level rows) from VEP, optionally weighted by `weight_col`."""
    # --> inputs <-- #
    # need a VEP executable
    if VEP_path is None:
        logger.error("VEP path not provided")
        return genes_df

    # pick the coordinate columns for the build
    if build.upper() in ["GRCH38", "HG38"]:
        chr_col = "chr_hg38_ensembl"
        pos_col = "pos_hg38"
        ref_col = "ref_hg38"
        alt_col = "alt_hg38"
        build_param = "GRCh38"
    else:
        chr_col = "chr_hg19_ensembl"
        pos_col = "pos_hg19"
        ref_col = "ref_hg19"
        alt_col = "alt_hg19"
        build_param = "GRCh37"

    # convert to pandas
    if isinstance(variants_df, pl.DataFrame):
        creds = variants_df.to_pandas()
    else:
        creds = variants_df.copy()

    if isinstance(genes_df, pl.DataFrame):
        genes_df = genes_df.to_pandas()

    # check the required columns
    required_cols = ["locus_id", chr_col, pos_col, ref_col, alt_col]
    for col in required_cols:
        if col not in creds.columns:
            logger.error(f"Required column {col} not found in input data")
            return genes_df

    # build a chr:pos:ref:alt variant id
    creds["variant_id"] = creds.apply(
        lambda row: f"{row[chr_col]}:{row[pos_col]}:{row[ref_col]}:{row[alt_col]}",
        axis=1,
    )

    # --> temp files <-- #
    # create the output dir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # unique per-run temp paths inside outdir
    tmpname = tempfile.mktemp(dir=outdir)
    tmp_file_path = f"{tmpname}.vcf"
    output_file_path = f"{tmpname}_output.vcf"
    variant_map_path = f"{tmpname}_variant_map.csv"

    # --> VEP input, both allele orientations <-- #
    output_lines = creds.copy()

    # sort by chromosome and position
    output_lines = output_lines.sort_values([chr_col, pos_col])

    # one map row per (locus, variant), a proxy shared by several loci gets one row each
    variant_map = []
    vcf_lines = []
    variant_map_flipped = []
    vcf_lines_flipped = []

    for idx, row in output_lines.iterrows():
        vcf_id = f"{row[chr_col]}_{row[pos_col]}_{row[ref_col]}/{row[alt_col]}"
        vcf_id_flipped = (
            f"{row[chr_col]}_{row[pos_col]}_{row[alt_col]}/{row[ref_col]}"
        )
        # end spans the longer allele
        end_pos = row[pos_col] + len(row[ref_col]) - 1
        end_pos_flipped = (
            row[pos_col] + len(row[alt_col]) - 1
        )

        # one input line per strand
        for strand in ["+", "-"]:
            vcf_line = f"{row[chr_col]} {row[pos_col]} {end_pos} {row[ref_col]}/{row[alt_col]} {strand}"
            vcf_line_flipped = f"{row[chr_col]} {row[pos_col]} {end_pos_flipped} {row[alt_col]}/{row[ref_col]} {strand}"
            vcf_lines.append(vcf_line)
            vcf_lines_flipped.append(vcf_line_flipped)

        # weight of this locus, 1.0 when the row carries none
        weight = 1.0
        if weight_col is not None and weight_col in row and not pd.isna(row[weight_col]):
            weight = row[weight_col]

        # map rows, flipped has ref and alt swapped
        variant_map.append({
            "vcf_id": vcf_id,
            "idx": idx,
            "variant_id": row["variant_id"],
            "locus_id": row["locus_id"],
            "chr": row[chr_col],
            "pos": row[pos_col],
            "ref": row[ref_col],
            "alt": row[alt_col],
            "scaling_factor": weight,
        })
        variant_map_flipped.append({
            "vcf_id": vcf_id_flipped,
            "idx": idx,
            "variant_id": row["variant_id"],
            "locus_id": row["locus_id"],
            "chr": row[chr_col],
            "pos": row[pos_col],
            "ref": row[alt_col],
            "alt": row[ref_col],
            "scaling_factor": weight,
        })
    # write the maps
    variant_map = pd.DataFrame(variant_map)
    variant_map_flipped = pd.DataFrame(variant_map_flipped)
    variant_map.to_csv(variant_map_path, index=False)
    variant_map_flipped.to_csv(variant_map_path + "_flipped", index=False)

    # dedupe in order and write the inputs
    vcf_lines = list(dict.fromkeys(vcf_lines))
    vcf_lines_flipped = list(dict.fromkeys(vcf_lines_flipped))
    with open(tmp_file_path, "w") as file:
        file.write("\n".join(vcf_lines))
    with open(tmp_file_path + "_flipped", "w") as file:
        file.write("\n".join(vcf_lines_flipped))

    # --> run VEP on both files <-- #
    try:
        run_vep_within_environment(
            input_file=tmp_file_path,
            output_file=output_file_path,
            vep_path=VEP_path,
            vep_fa_cache=VEP_fa_cache,
            vep_dir_cache=VEP_dir_cache,
            build=build_param,
            vep_env=VEP_env,
            env_manager=VEP_env_manager,
            fork=VEP_fork,
        )
        logger.info("- - VEP completed successfully for unflipped variants")
        run_vep_within_environment(
            input_file=f"{tmp_file_path}_flipped",
            output_file=f"{output_file_path}_flipped",
            vep_path=VEP_path,
            vep_fa_cache=VEP_fa_cache,
            vep_dir_cache=VEP_dir_cache,
            build=build_param,
            vep_env=VEP_env,
            env_manager=VEP_env_manager,
            fork=VEP_fork,
        )
        logger.info("- - VEP completed successfully for flipped variants")

    except Exception as e:
        logger.error(f"Error running VEP: {e}")
        return genes_df

    # --> read VEP output <-- #
    # skip the ## meta lines
    data = []
    data_flipped = []
    try:
        with open(output_file_path, "r") as file:
            for line in file:
                if not line.startswith("##"):
                    data.append(line.strip().split("\t"))
        with open(output_file_path + "_flipped", "r") as file:
            for line in file:
                if not line.startswith("##"):
                    data_flipped.append(line.strip().split("\t"))
    except Exception as e:
        logger.error(f"Error reading VEP output: {e}")
        return genes_df

    # stop when both files hold a header at most
    if (
        len(data) <= 1 and len(data_flipped) <= 1
    ):
        logger.warning("No data found in VEP output")
        return genes_df

    # --> parse VEP output <-- #
    try:
        # first row is the header
        vep_df = (
            pd.DataFrame(data[1:], columns=data[0]) if len(data) > 1 else pd.DataFrame()
        )
        vep_df_flipped = (
            pd.DataFrame(data_flipped[1:], columns=data_flipped[0])
            if len(data_flipped) > 1
            else pd.DataFrame()
        )

        # --> unflipped variants <-- #
        if not vep_df.empty:
            # IMPACT from the Extra field
            vep_df["impact"] = vep_df["Extra"].str.extract(r"IMPACT=(\w+);")
            vep_df["impact"] = vep_df["impact"].fillna("None")

            # CADD scores from the Extra field
            vep_df["cadd_phred"] = vep_df["Extra"].str.extract(r"CADD_PHRED=([^;]+)")
            vep_df["cadd_raw"] = vep_df["Extra"].str.extract(r"CADD_RAW=([^;]+)")

            # cast CADD to float, NaN when missing
            vep_df["cadd_phred"] = pd.to_numeric(vep_df["cadd_phred"], errors="coerce")
            vep_df["cadd_raw"] = pd.to_numeric(vep_df["cadd_raw"], errors="coerce")

            # impact to score
            VEP_dict = {
                "HIGH": 1,
                "MODERATE": 0.6,
                "LOW": 0.4,
                "MODIFIER": 0.1,
                "None": 0,
            }
            vep_df["consequence_score"] = vep_df["impact"].apply(lambda x: VEP_dict[x])

            # map back to every (locus, variant) the uploaded variation belongs to
            vep_df = _map_to_loci(vep_df, variant_map, keep_weight=bool(use_weights and weight_col))

            # Apply r2 weighting if requested #TODO something isn t right here
            if use_weights and weight_col:
                vep_df["consequence_score"] = vep_df["consequence_score"] * vep_df["scaling_factor"]

            # rename the Ensembl gene column
            vep_df = vep_df.rename(columns={"Gene": "ensembl_gene_id"})

        # --> flipped variants, same steps <-- #
        if not vep_df_flipped.empty:
            # IMPACT from the Extra field
            vep_df_flipped["impact"] = vep_df_flipped["Extra"].str.extract(
                r"IMPACT=(\w+);"
            )
            vep_df_flipped["impact"] = vep_df_flipped["impact"].fillna("None")

            # CADD scores from the Extra field
            vep_df_flipped["cadd_phred"] = vep_df_flipped["Extra"].str.extract(
                r"CADD_PHRED=([^;]+)"
            )
            vep_df_flipped["cadd_raw"] = vep_df_flipped["Extra"].str.extract(
                r"CADD_RAW=([^;]+)"
            )

            # cast CADD to float, NaN when missing
            vep_df_flipped["cadd_phred"] = pd.to_numeric(
                vep_df_flipped["cadd_phred"], errors="coerce"
            )
            vep_df_flipped["cadd_raw"] = pd.to_numeric(
                vep_df_flipped["cadd_raw"], errors="coerce"
            )

            # impact to score
            vep_df_flipped["consequence_score"] = vep_df_flipped["impact"].apply(
                lambda x: VEP_dict[x]
            )

            # map back through the flipped map
            vep_df_flipped = _map_to_loci(
                vep_df_flipped, variant_map_flipped, keep_weight=bool(use_weights and weight_col)
            )

            # weight the consequence score
            if use_weights and weight_col:
                vep_df_flipped["consequence_score"] = (
                    vep_df_flipped["consequence_score"] * vep_df_flipped["scaling_factor"]
                )

            # rename the Ensembl gene column
            vep_df_flipped = vep_df_flipped.rename(columns={"Gene": "ensembl_gene_id"})

        # --> combine and summarise per gene <-- #
        # stack whichever orientations returned rows
        if not vep_df.empty and not vep_df_flipped.empty:
            combined_vep_df = pd.concat([vep_df, vep_df_flipped], ignore_index=True)
        elif not vep_df.empty:
            combined_vep_df = vep_df
        elif not vep_df_flipped.empty:
            combined_vep_df = vep_df_flipped
        else:
            logger.warning("No valid data found in either VEP output")
            return genes_df
        
        # drop rows both orientations returned
        if {"locus_id", "variant_id", "ensembl_gene_id"}.issubset(set(combined_vep_df.columns)):
            combined_vep_df = combined_vep_df.drop_duplicates(
                subset=["locus_id", "variant_id", "ensembl_gene_id"], keep="first"
            )
        
        # gene-level summaries
        agg_dict = {
            "consequence_score": ["max", "sum",],
            "cadd_phred": ["max", "mean","sum"],
            "cadd_raw": ["max", "mean","sum"],
        }
        if use_weights and "scaling_factor" in combined_vep_df.columns:
            agg_dict["scaling_factor"] = ["sum"]        


        # aggregate per (locus, gene)
        gene_summary = combined_vep_df.groupby(["locus_id", "ensembl_gene_id"]).agg(agg_dict)

        # flatten the MultiIndex into the output names
        new_cols = []
        for col, stat in gene_summary.columns:
            if col == "consequence_score" and stat == "max":
                new_cols.append("VEP_max")
            elif col == "consequence_score" and stat == "sum":
                new_cols.append("VEP_sum")
            elif col == "cadd_phred" and stat == "max":
                new_cols.append("CADD_phred_max")
            elif col == "cadd_phred" and stat == "mean":
                new_cols.append("CADD_phred_mean")
            elif col == "cadd_phred" and stat == "sum":
                new_cols.append("CADD_phred_sum")
            elif col == "cadd_raw" and stat == "max":
                new_cols.append("CADD_raw_max")
            elif col == "cadd_raw" and stat == "mean":
                new_cols.append("CADD_raw_mean")
            elif col == "cadd_raw" and stat == "sum":
                new_cols.append("CADD_raw_sum")
            elif col == "scaling_factor" and stat == "sum":
                new_cols.append("weight_sum")
            else:
                new_cols.append(f"{col}_{stat}")
                
        gene_summary.columns = new_cols
        gene_summary = gene_summary.reset_index()
        # weighted means when weights are available
        # `consequence_score` is already weighted, so VEP_sum / weight_sum is the weighted mean
        if use_weights and "weight_sum" in gene_summary.columns:
            gene_summary["VEP_weighted_mean"] = gene_summary["VEP_sum"] / gene_summary["weight_sum"]
            gene_summary["CADD_phred_weighted_mean"] = gene_summary["CADD_phred_sum"] / gene_summary["weight_sum"]
            gene_summary["CADD_raw_weighted_mean"] = gene_summary["CADD_raw_sum"] / gene_summary["weight_sum"]

        # back to polars when the input was polars
        if isinstance(variants_df, pl.DataFrame):
            combined_vep_df = pl.from_pandas(combined_vep_df)
            gene_summary = pl.from_pandas(gene_summary)

        return gene_summary, combined_vep_df

    except Exception as e:
        logger.error(f"Error processing VEP results: {e}")
        return genes_df

    finally:
        # remove the temp files
        for f in [
            tmp_file_path,
            tmp_file_path + "_flipped",
            output_file_path,
            output_file_path + "_flipped",
            output_file_path + "_summary.html",
            output_file_path + "_flipped_summary.html",
            output_file_path + "_warnings.txt",
            output_file_path + "_flipped_warnings.txt",
            variant_map_path,
            variant_map_path + "_flipped",
        ]:
            if os.path.exists(f):
                os.remove(f)


################################################
####                vep run                 ####
################################################


def _activation_prefix(env_manager, env):
    """Shell command that activates `env` with `env_manager`."""
    if env_manager not in ENV_MANAGERS:
        raise ValueError(
            f"Unknown vep.env_manager {env_manager!r}, "
            f"expected one of {sorted(ENV_MANAGERS)}."
        )
    return ENV_MANAGERS[env_manager].format(env=env)

def run_vep_within_environment(
    input_file,
    output_file,
    vep_path,
    vep_fa_cache,
    vep_dir_cache,
    build,
    vep_env=None,
    env_manager="micromamba",
    fork=8,
):
    """Run VEP offline with the CADD and SpliceRegion plugins inside `vep_env`."""
    # need an env to activate
    if vep_env is None:
        raise ValueError(
            "No conda environment given for VEP. Set vep.env "
            "in your config (see configs/base_config.yaml)."
        )
    # normalise the assembly name
    if build.upper() == "GRCH38" or build.upper() == "HG38":
        build = "GRCh38"
    else:
        logger.info("Only supporting GRCh38 build for VEP at this time")

    # TODO fix this
    cadd_path = os.path.join(
        vep_dir_cache,
        "vep_plugins/CADD/krishna.gs.washington.edu/download/CADD/v1.6/GRCh38",
    )
    dir_plugis = os.path.join(vep_dir_cache, ".vep/Plugins38")
    vep_command = f"{vep_path} -i {input_file} -o {output_file} --fasta {vep_fa_cache} --dir_cache {vep_dir_cache} --assembly {build}  --dir_plugins {dir_plugis} --plugin CADD,{cadd_path}/whole_genome_SNVs.tsv.gz,{cadd_path}/gnomad.genomes.r3.0.indel.tsv.gz --plugin SpliceRegion --offline --fork {fork} --force_overwrite"
    # prefix with the env activation
    full_command = f"{_activation_prefix(env_manager, vep_env)} && {vep_command}"

    # run under bash
    subprocess.run(full_command, shell=True, check=True, executable="/bin/bash")
    return


################################################
####              entry point               ####
################################################


def get_VEP_annotations(
    variant_data,
    genes_in_region,
    config,
    use_scaling_factor=True,
):
    """VEP (gene summary, variant-level rows) for all locus variants, weighted by `scaling_factor`."""
    logger.info("Starting VEP annotation using command line VEP")

    # weight column
    weight_col = "scaling_factor"

    # run command-line VEP with the config settings
    return process_variants_with_cmd_vep(
        variants_df=variant_data,
        genes_df=genes_in_region,
        weight_col=weight_col if use_scaling_factor else None,
        use_weights=use_scaling_factor,
        build=config.parameters.genome_build,
        VEP_path=config.vep.executable,
        VEP_fa_cache=config.vep.fasta_file,
        VEP_dir_cache=config.vep.cache_dir,
        VEP_env=config.vep.env,
        VEP_env_manager=config.vep.env_manager,
        VEP_fork=config.vep.fork,
        # temp dir follows the run cache, so parallel jobs never share one
        outdir=str(get_cache_dir() / "vep_tempdir"),
    )
