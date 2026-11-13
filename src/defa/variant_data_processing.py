################################################
#### Variant data: leads, proxies, liftover ####
################################################
"""Load leads and proxies, lift over to the missing build and build marker names in both builds."""

import polars as pl
from loguru import logger
import sys
from pathlib import Path
from pyliftover import LiftOver

from .input_schema import read_input_table

sys.path.append(str(Path(__file__).parent.parent.parent / "helper_scripts"))


################################################
####      coordinates and marker names      ####
################################################


def add_liftover_positions(
    df: pl.DataFrame, from_build: str, to_build: str
) -> pl.DataFrame:
    """Lift `pos_<from_build>` to `pos_<to_build>` and copy chr, ref and alt across."""
    chr_col = f"chr_{from_build}"
    pos_col = f"pos_{from_build}"
    target_pos_col = f"pos_{to_build}"

    # skip when the source columns are missing
    if chr_col not in df.columns or pos_col not in df.columns:
        logger.warning(f"Columns {chr_col} or {pos_col} missing. Skipping liftover.")
        return df

    logger.info(f"Lifting over from {from_build} to {to_build}...")

    # set up LiftOver
    # pyliftover derives the chain filename itself and downloads it on first use
    try:
        lo = LiftOver(from_build, to_build)
    except Exception as e:
        logger.error(f"Failed to initialize LiftOver {from_build}->{to_build}: {e}")
        return df

    # --> per-row mapping <-- #
    def perform_liftover_map(coords: dict) -> int | None:
        """Lifted position for one (chr, pos) struct, None when it does not map."""
        chrom = coords.get(chr_col)
        pos = coords.get(pos_col)

        if chrom is None or pos is None:
            return None

        try:
            pos = int(pos)
        except (ValueError, TypeError):
            return None

        # normalise to the chr-prefixed UCSC name
        chrom_str = str(chrom)
        if not chrom_str.startswith("chr"):
            if chrom_str in ["MT", "M"]:
                chrom_str = "chrM"
            else:
                chrom_str = "chr" + chrom_str

        # take the first mapped position
        try:
            result = lo.convert_coordinate(chrom_str, pos)
            if result:
                return result[0][1]
        except Exception:
            return None
        return None

    # --> lift and copy <-- #
    df_lifted = df.with_columns(
        [
            pl.struct([chr_col, pos_col])
            .map_elements(
                perform_liftover_map,
                return_dtype=pl.Int64,
                skip_nulls=False,
            )
            .alias(target_pos_col),
            # chr, ref and alt carry over unchanged
            pl.col(f"chr_{from_build}").alias(f"chr_{to_build}"),
            pl.col(f"ref_{from_build}").alias(f"ref_{to_build}"),
            pl.col(f"alt_{from_build}").alias(f"alt_{to_build}"),
        ]
    )

    return df_lifted


def create_marker_name(
    variant_df: pl.DataFrame,
    chr_column_name: str = "CHR",
    pos_column_name: str = "POS",
    allele_0_column_name: str = "REF",
    allele_1_column_name: str = "ALT",
    desired_output_format: str = "ucsc",
    genome_build: str = "hg38",
    check_nans=False,
) -> pl.Series:
    """Marker ids as chr:pos_A1_A2, alleles uppercased and sorted alphabetically."""
    # optionally fail on any null
    if check_nans:
        assert variant_df.null_count().pipe(sum).item() == 0, (
            "Variant Dataset Contain NaNs, please check!"
        )

    # select and cast the four columns
    marker_name_data = variant_df.select(
        [
            pl.col(chr_column_name).cast(pl.String),
            pl.col(pos_column_name).cast(pl.Int64),
            pl.col(allele_0_column_name).cast(pl.String),
            pl.col(allele_1_column_name).cast(pl.String),
        ]
    )

    # replace 23/24 with X/Y
    marker_name_data = marker_name_data.with_columns(
        [
            pl.col(chr_column_name)
            .str.replace("23", "X")
            .str.replace("24", "Y")
            .alias(chr_column_name)
        ]
    )

    # --> UCSC and Ensembl chromosome names <-- #
    # add the chr prefix when it is missing, strip it for the Ensembl name otherwise
    if (
        not marker_name_data.select(pl.col(chr_column_name).str.starts_with("chr"))
        .to_series()
        .all()
    ):
        marker_name_data = marker_name_data.with_columns(
            [
                (pl.lit("chr") + pl.col(chr_column_name)).alias("ucsc_chr_name"),
                pl.col(chr_column_name).alias("ensembl_chr_name"),
            ]
        )
    else:
        marker_name_data = marker_name_data.with_columns(
            pl.col(chr_column_name).alias("ucsc_chr_name"),
            pl.when(pl.col(chr_column_name).str.contains("chrX"))
            .then(pl.lit("X"))
            .when(pl.col(chr_column_name).str.contains("chrY"))
            .then(pl.lit("Y"))
            .otherwise(pl.col(chr_column_name).str.extract(r"(\d+)", 1))
            .alias("ensembl_chr_name"),
        )

    # pick the chromosome naming for the output
    chromosme_type = (
        "ucsc_chr_name" if desired_output_format == "ucsc" else "ensembl_chr_name"
    )
    marker_col_name = f"marker_name_{genome_build}"

    # --> marker ids <-- #
    marker_name_data = marker_name_data.with_columns(
        [
            pl.concat_str([chromosme_type, pos_column_name], separator=":").alias(
                marker_col_name
            ),
            # order alleles alphabetically
            pl.concat_str(
                [
                    pl.when(
                        pl.col(allele_0_column_name).str.to_uppercase()
                        < pl.col(allele_1_column_name).str.to_uppercase()
                    )
                    .then(pl.col(allele_0_column_name).str.to_uppercase())
                    .otherwise(pl.col(allele_1_column_name).str.to_uppercase()),
                    pl.when(
                        pl.col(allele_0_column_name).str.to_uppercase()
                        < pl.col(allele_1_column_name).str.to_uppercase()
                    )
                    .then(pl.col(allele_1_column_name).str.to_uppercase())
                    .otherwise(pl.col(allele_0_column_name).str.to_uppercase()),
                ],
                separator="_",
            ).alias("ordered_alleles"),
        ]
    )

    # append the ordered alleles
    marker_name_data = marker_name_data.with_columns(
        pl.concat_str([marker_col_name, "ordered_alleles"], separator="_").alias(
            marker_col_name
        )
    )

    return marker_name_data[marker_col_name]


def add_ensembl_columns(df: pl.DataFrame, build: str) -> pl.DataFrame:
    """Add `chr_<build>_ensembl`, the chromosome without the chr prefix and M as MT."""
    chr_col = f"chr_{build}"
    ens_col = f"chr_{build}_ensembl"

    # map X, Y and M, keep the digits otherwise
    if chr_col in df.columns:
        return df.with_columns(
            pl.when(pl.col(chr_col).str.contains("(?i)chrX"))
            .then(pl.lit("X"))
            .when(pl.col(chr_col).str.contains("(?i)chrY"))
            .then(pl.lit("Y"))
            .when(pl.col(chr_col).str.contains("(?i)chrM"))
            .then(pl.lit("MT"))
            .otherwise(pl.col(chr_col).str.extract(r"(\d+)", 1))
            .alias(ens_col)
        )
    return df


################################################
####           variant processing           ####
################################################


def process_variants_dataframe(
    df: pl.DataFrame, scaling_factor_threshold: float = 0.0
) -> pl.DataFrame:
    """Filter by `scaling_factor`, lift over to the missing build and add marker names in both builds."""

    # --> scaling factor <-- #
    # rename r2 when scaling_factor is absent
    if "scaling_factor" not in df.columns:
        if "r2" in df.columns:
            logger.info("Column 'r2' detected. Renaming to 'scaling_factor'.")
            df = df.rename({"r2": "scaling_factor"})
        else:
            logger.warning("No score column found. 'scaling_factor' will be missing.")

    # filter by threshold
    if "scaling_factor" in df.columns:
        df = df.with_columns(pl.col("scaling_factor").cast(pl.Float64))
        if scaling_factor_threshold > 0:
            df = df.filter(pl.col("scaling_factor") >= scaling_factor_threshold)

    # --> source and destination build <-- #
    has_hg38 = "pos_hg38" in df.columns
    has_hg19 = "pos_hg19" in df.columns

    if not has_hg38 and not has_hg19:
        raise ValueError("DataFrame must contain either pos_hg38 or pos_hg19 columns.")

    # hg38 is the source whenever it is present
    src_build = "hg38" if has_hg38 else "hg19"
    dst_build = "hg19" if src_build == "hg38" else "hg38"

    logger.info(f"Processing variants. Source: {src_build} -> Target: {dst_build}")

    # --> source build <-- #
    # cast position, add start and end
    df = df.with_columns(
        [
            pl.col(f"pos_{src_build}").cast(pl.Int64),
            pl.col(f"pos_{src_build}").alias(f"start_{src_build}"),
            pl.col(f"pos_{src_build}").alias(f"end_{src_build}"),
        ]
    )

    df = add_ensembl_columns(df, src_build)

    # source marker name
    df = df.with_columns(
        create_marker_name(
            df,
            chr_column_name=f"chr_{src_build}",
            pos_column_name=f"pos_{src_build}",
            allele_0_column_name=f"ref_{src_build}",
            allele_1_column_name=f"alt_{src_build}",
            genome_build=src_build,
            desired_output_format="ucsc",
        ).alias(f"marker_name_{src_build}")
    )

    # --> destination build <-- #
    if not (has_hg38 and has_hg19):
        # lift over when the destination is missing
        df = add_liftover_positions(df, from_build=src_build, to_build=dst_build)

        # drop unmapped variants, then add start, end and marker name
        if f"pos_{dst_build}" in df.columns:
            df = df.drop_nulls(subset=[f"pos_{dst_build}"])
            df = add_ensembl_columns(df, dst_build)
            df = df.with_columns(
                [
                    pl.col(f"pos_{dst_build}").alias(f"start_{dst_build}"),
                    pl.col(f"pos_{dst_build}").alias(f"end_{dst_build}"),
                ]
            )
            df = df.with_columns(
                create_marker_name(
                    df,
                    chr_column_name=f"chr_{dst_build}",
                    pos_column_name=f"pos_{dst_build}",
                    allele_0_column_name=f"ref_{dst_build}",
                    allele_1_column_name=f"alt_{dst_build}",
                    genome_build=dst_build,
                    desired_output_format="ucsc",
                ).alias(f"marker_name_{dst_build}")
            )
        else:
            logger.error(f"Liftover failed. {dst_build} columns will be missing.")
    else:
        # both builds given, fill in any missing helper columns
        for b in ["hg19", "hg38"]:
            df = add_ensembl_columns(df, b)
            if f"marker_name_{b}" not in df.columns:
                df = df.with_columns(
                    create_marker_name(
                        df,
                        chr_column_name=f"chr_{b}",
                        pos_column_name=f"pos_{b}",
                        allele_0_column_name=f"ref_{b}",
                        allele_1_column_name=f"alt_{b}",
                        genome_build=b,
                        desired_output_format="ucsc",
                    ).alias(f"marker_name_{b}")
                )
            if f"start_{b}" not in df.columns:
                df = df.with_columns(
                    [
                        pl.col(f"pos_{b}").alias(f"start_{b}"),
                        pl.col(f"pos_{b}").alias(f"end_{b}"),
                    ]
                )

    # --> output columns <-- #
    # keep the ones present
    cols_to_keep = [
        "locus_id",
        "variant_id",
        "marker_name_hg38",
        "marker_name_hg19",
        "scaling_factor",
        # hg38 columns
        "chr_hg38",
        "pos_hg38",
        "start_hg38",
        "end_hg38",
        "ref_hg38",
        "alt_hg38",
        "chr_hg38_ensembl",
        # hg19 columns
        "chr_hg19",
        "pos_hg19",
        "start_hg19",
        "end_hg19",
        "ref_hg19",
        "alt_hg19",
        "chr_hg19_ensembl",
    ]
    final_cols = [c for c in cols_to_keep if c in df.columns]

    return df.select(final_cols)


################################################
####                loading                 ####
################################################


def validate_input_file(df: pl.DataFrame) -> pl.DataFrame:
    """Rename r2 to `scaling_factor` and raise ValueError on any missing id or coordinate column."""

    # --> score column <-- #
    if "scaling_factor" in df.columns:
        pass
    elif "r2" in df.columns:
        df = df.rename({"r2": "scaling_factor"})
    else:
        raise ValueError("Missing Score Column: Input must contain 'r2'.")

    # --> id columns <-- #
    missing_ids = [
        col for col in ["locus_id", "variant_id"] if col not in df.columns
    ]
    if missing_ids:
        raise ValueError(f"Missing ID Columns: Input must contain {missing_ids}")

    # --> coordinate columns, hg38 OR hg19 <-- #
    has_hg38 = "pos_hg38" in df.columns
    has_hg19 = "pos_hg19" in df.columns

    if not (has_hg38 or has_hg19):
        raise ValueError(
            "Missing Coordinates: Input must contain either 'pos_hg38' or 'pos_hg19'."
        )

    # chr, ref and alt for whichever build is present
    missing_genome_cols = []

    if has_hg38:
        required = ["chr_hg38", "ref_hg38", "alt_hg38"]
        missing = [c for c in required if c not in df.columns]
        if missing:
            missing_genome_cols.extend(missing)

    if has_hg19:
        required = ["chr_hg19", "ref_hg19", "alt_hg19"]
        missing = [c for c in required if c not in df.columns]
        if missing:
            missing_genome_cols.extend(missing)

    if missing_genome_cols:
        raise ValueError(
            f"Missing Genomic Data: detected position column but missing {missing_genome_cols}"
        )

    return df


# coordinate columns a lead or proxy file may carry
GENOME_COLUMNS = [
    "chr_hg38", "pos_hg38", "ref_hg38", "alt_hg38",
    "chr_hg19", "pos_hg19", "ref_hg19", "alt_hg19",
]


def load_lead_and_proxy(
    lead_path: str, proxy_path: str, scaling_factor_threshold: float = 0.0
) -> pl.DataFrame:
    """Stack the lead and proxy files into one row per (locus, variant) and process it."""
    logger.info(f"Loading lead file {lead_path} and proxy file {proxy_path}")
    lead_df  = read_input_table(lead_path, "lead").drop_nulls()
    proxy_df = read_input_table(proxy_path, "proxy").drop_nulls()

    # a lead is its own variant at full weight
    lead_df = validate_input_file(lead_df.with_columns(
        pl.col("locus_id").alias("variant_id"),
        pl.lit(1.0).alias("scaling_factor"),
    ))
    proxy_df = validate_input_file(proxy_df.rename({"proxy_id": "variant_id"}))

    # keep the ids, the score, and the coordinate columns both files share
    shared = [c for c in GENOME_COLUMNS if c in lead_df.columns and c in proxy_df.columns]
    keep   = ["locus_id", "variant_id", "scaling_factor"] + shared
    combined_df = pl.concat(
        [lead_df.select(keep), proxy_df.select(keep)], how="vertical_relaxed"
    )
    # a variant listed twice for one locus keeps its higher score, so a lead keeps 1.0
    combined_df = combined_df.sort("scaling_factor", descending=True).unique(
        subset=["locus_id", "variant_id"], keep="first"
    )
    return process_variants_dataframe(combined_df, scaling_factor_threshold)


def collect_variant_data(config) -> pl.DataFrame:
    """Variant data from the lead and proxy files together."""
    dataset   = config.dataset
    threshold = config.parameters.scaling_factor_threshold

    if not (dataset.lead_variant_file and dataset.proxy_file):
        raise ValueError(
            "no variant input configured: set dataset.lead_variant_file AND dataset.proxy_file"
        )
    variant_data = load_lead_and_proxy(
        dataset.lead_variant_file, dataset.proxy_file, threshold
    )

    logger.info(f"Collected {variant_data.height} variants")
    return variant_data
