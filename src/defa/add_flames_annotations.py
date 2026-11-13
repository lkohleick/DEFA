################################################
####   FLAMES: regulatory and QTL evidence  ####
################################################
"""FLAMES annotations per (locus, gene): ABC, HACER, QTL colocalisation, Hi-C and enhancer-gene links."""

import os

import polars as pl
from loguru import logger
from tqdm import tqdm


# the 17 Javierre cell types, `sum` and `max` are the horizontal aggregates of these
JAVIERRE_CELL_TYPES = [
    "Mon", "Mac0", "Mac1", "Mac2", "Neu", "MK", "EP", "Ery", "FoeT",
    "nCD4", "tCD4", "aCD4", "naCD4", "nCD8", "tCD8", "nB", "tB",
]


################################################
####       per gene annotation files        ####
################################################


def _list_gene_files(directory):
    """Map every `ENSG<id>.parquet` name under `directory` to its path."""
    paths = {}
    for entry in os.listdir(directory):
        full = os.path.join(directory, entry)

        # the flat layout, which every dataset but eQTLGen uses
        if entry.startswith("ENSG") and entry.endswith(".parquet"):
            paths[entry] = full
            continue

        # eQTLGen shards its per-gene files into chr1 ... chr22
        if entry.startswith("chr") and os.path.isdir(full):
            for name in os.listdir(full):
                if name.startswith("ENSG") and name.endswith(".parquet"):
                    paths[name] = os.path.join(full, name)

    return paths


def _read_gene_parquet_files(directory, gene_ids, gene_col=None):
    """Read the per-gene parquet files under `directory` for the genes that have one."""
    # list the directory once, a stat per gene is thousands of round trips on a network mount
    available = _list_gene_files(directory)

    # an empty directory means the config points somewhere wrong
    if not available:
        raise FileNotFoundError(f"no ENSG<id>.parquet files found in {directory}")

    wanted = sorted({f"{gene}.parquet" for gene in gene_ids} & set(available))
    logger.info(f"Loading gene files from {directory}")

    if not wanted:
        probe = available[min(available)]
        schema = dict(pl.scan_parquet(probe).collect_schema())
        # eQTLGen keeps the gene only in the file name, so the column has to be added back
        if gene_col:
            schema[gene_col] = pl.String
        return pl.DataFrame(schema=schema)

    frames = []
    for name in wanted:
        frame = pl.scan_parquet(available[name])
        if gene_col:
            frame = frame.with_columns(pl.lit(name[:-8]).alias(gene_col))
        frames.append(frame)

    # `vertical_relaxed` coerces the int64 versus float64 drift between per gene files
    return pl.concat(frames, how="vertical_relaxed").collect()


def _marker_name_hg38(source):
    """Build `chr:pos_A1_A2` with the alleles alphabetical, as the combined files store it."""
    # the source field looks like chr20_50049160_G_A_b38 or chrX_100649875_A_G
    parts = source.str.replace("_b38$", "").str.split("_")
    allele_1 = parts.list.get(2)
    allele_2 = parts.list.get(3)
    return pl.concat_str(
        [
            parts.list.get(0),
            pl.lit(":"),
            parts.list.get(1),
            pl.lit("_"),
            pl.min_horizontal(allele_1, allele_2),
            pl.lit("_"),
            pl.max_horizontal(allele_1, allele_2),
        ]
    ).alias("marker_name_hg38")


################################################
####             abc and hacer              ####
################################################


def get_ABC_EP_annotations(variant_data, genes_in_region, config, prob_col=None):
    """ABC enhancer-promoter scores of the variants inside each candidate gene's elements, sum, max and mean per (locus, gene)."""
    dataset_name = "ABC_EP"
    gene_col_dataset = "ENSG"
    score_col_dataset = "ABC.Score"
    # load ABC EP data from the per gene files
    abc_ep_data = _read_gene_parquet_files(
        config.paths.gene_regulation.abc_ep_data,
        genes_in_region["ensembl_gene_id"],
    )
    abc_ep_data = abc_ep_data.rename({gene_col_dataset: "ensembl_gene_id"})
    if not prob_col:
        prob_col = "prob_col"
        variant_data = variant_data.with_columns(
            pl.lit(1).alias(prob_col).cast(pl.Float64)
        )
    # log the load
    logger.info(f"Loaded ABC EP data: {abc_ep_data.shape}")

    # split the variants by locus
    locus_groups_pl = variant_data.partition_by("locus_id")
    # annotate each locus
    loci_annotations = []

    for locus_var_df in tqdm(locus_groups_pl, desc="Processing loci"):
        list_of_genes_in_region = genes_in_region.filter(
            pl.col("locus_id").is_in(locus_var_df["locus_id"].to_list())
        )["ensembl_gene_id"]
        annotations_genes = abc_ep_data.filter(
            pl.col("ensembl_gene_id").is_in(list_of_genes_in_region)
        )
        if annotations_genes.shape[0] == 0:
            continue
        # do a cross join between the subsetted annotations and the locus_id information
        cross_df = annotations_genes.select(
            [
                pl.col("chr"),
                pl.col("start"),
                pl.col("end"),
                pl.col("ensembl_gene_id"),
                pl.col(score_col_dataset),
            ]
        ).join(
            locus_var_df.select(
                [
                    pl.col("locus_id"),
                    pl.col("chr_hg19"),
                    pl.col("pos_hg19"),
                    pl.col(prob_col),
                ]
            ),
            how="cross",
        )
        # filter for exact matches in region
        matches_df = cross_df.filter(
            (pl.col("chr") == pl.col("chr_hg19"))
            & (pl.col("pos_hg19") >= pl.col("start"))
            & (pl.col("pos_hg19") <= pl.col("end"))
        )
        # calculate the scaled measurement (with probability or r2)
        matches_df = matches_df.with_columns(
            (pl.col(score_col_dataset) * pl.col(prob_col)).alias("scaled_measurement")
        )
        # keep keys and the measurement
        matches_df = matches_df.select(
            [
                pl.col("locus_id"),
                pl.col("ensembl_gene_id"),
                pl.col("scaled_measurement"),
            ]
        )
        # aggregate the data with sum, max, mean
        results_df = matches_df.group_by(["locus_id", "ensembl_gene_id"]).agg(
            pl.sum("scaled_measurement").alias(f"{dataset_name}_sum"),
            pl.max("scaled_measurement").alias(f"{dataset_name}_max"),
            pl.mean("scaled_measurement").alias(f"{dataset_name}_mean"),
        )

        loci_annotations.append(results_df)
    # stack the loci
    results_df = pl.concat(loci_annotations)
    return results_df


def get_ABC_CRISPR_annotations(variant_data, genes_in_region, config, prob_col=None):
    """ABC CRISPR-validated scores per (locus, gene), with sum, max, mean and the weighted mean."""
    dataset_name = "ABC_CRISPR"
    gene_col_dataset = "ENSG"
    score_col_dataset = "ABC_max"
    # load ABC CRISPR data from the per gene files
    abc_ep_data = _read_gene_parquet_files(
        config.paths.gene_regulation.abc_crispr_data,
        genes_in_region["ensembl_gene_id"],
    )
    abc_ep_data = abc_ep_data.rename({gene_col_dataset: "ensembl_gene_id"})
    if not prob_col:
        prob_col = "prob_col"
        variant_data = variant_data.with_columns(
            pl.lit(1).alias(prob_col).cast(pl.Float64)
        )
    # log the load
    logger.info(f"Loaded ABC CRISPR data: {abc_ep_data.shape}")

    # split the variants by locus
    locus_groups_pl = variant_data.partition_by("locus_id")
    # annotate each locus
    loci_annotations = []

    for locus_var_df in tqdm(locus_groups_pl, desc="Processing loci"):
        list_of_genes_in_region = genes_in_region.filter(
            pl.col("locus_id") == locus_var_df["locus_id"][0]
        )["ensembl_gene_id"]
        annotations_genes = abc_ep_data.filter(
            pl.col("ensembl_gene_id").is_in(list_of_genes_in_region)
        )
        if annotations_genes.shape[0] == 0:
            continue
        # do a cross join between the subsetted annotations and the locus_id information
        cross_df = annotations_genes.select(
            [
                pl.col("chr"),
                pl.col("start"),
                pl.col("end"),
                pl.col("ensembl_gene_id"),
                pl.col(score_col_dataset),
            ]
        ).join(
            locus_var_df.select(
                [
                    pl.col("locus_id"),
                    pl.col("chr_hg19"),
                    pl.col("pos_hg19"),
                    pl.col(prob_col),
                ]
            ),
            how="cross",
        )
        # filter for exact matches in region
        matches_df = cross_df.filter(
            (pl.col("chr") == pl.col("chr_hg19"))
            & (pl.col("pos_hg19") >= pl.col("start"))
            & (pl.col("pos_hg19") <= pl.col("end"))
        )
        # calculate the scaled measurement (with probability or r2)
        matches_df = matches_df.with_columns(
            (pl.col(score_col_dataset) * pl.col(prob_col)).alias("scaled_measurement")
        )
        # keep keys and the measurement
        matches_df = matches_df.select(
            [
                pl.col("locus_id"),
                pl.col("ensembl_gene_id"),
                pl.col("scaled_measurement"),
                pl.col(prob_col)
            ]
        )
        # aggregate the data with sum, max, mean
        results_df = matches_df.group_by(["locus_id", "ensembl_gene_id"]).agg(
            pl.sum("scaled_measurement").alias(f"{dataset_name}_sum"),
            pl.max("scaled_measurement").alias(f"{dataset_name}_max"),
            pl.mean("scaled_measurement").alias(f"{dataset_name}_mean"),
            pl.col(prob_col).sum().alias("weight_sum"),
        )
        # proper weighted mean: (sum score * weight) / sum weight
        results_df = results_df.with_columns(
            (pl.col(f"{dataset_name}_sum") / pl.col("weight_sum")).alias(
                f"{dataset_name}_weighted_mean"
            )
        )
        loci_annotations.append(results_df)
    # stack the loci
    results_df = pl.concat(loci_annotations)
    return results_df


def get_HACER_annotations(variant_data, genes_in_region, config, prob_col=None):
    """HACER CAGE and PRO-seq/GRO-seq enhancer densities per (locus, gene), pooled over the F5, closest and 4D links."""
    logger.info(" - Starting HACER annotations")

    gene_col_dataset = "ensg_id"
    score_col_dataset = "density"
    if not prob_col:
        prob_col = "prob_col"
        variant_data = variant_data.with_columns(
            pl.lit(1).alias(prob_col).cast(pl.Float64)
        )
    hacer_results = {}
    # split the variants by locus
    locus_groups_pl = variant_data.partition_by("locus_id")

    def get_subhacer_annotations(list_of_locus_dfs, dataset, genes_in_region):
        """Weighted HACER densities of the variants inside each candidate gene's enhancers, one row per match."""
        loci_annotations = []
        for locus_var_df in tqdm(list_of_locus_dfs, desc="Processing loci"):
            list_of_genes_in_region = genes_in_region.filter(
                pl.col("locus_id") == locus_var_df["locus_id"][0]
            )["ensembl_gene_id"]
            annotations_genes = dataset.filter(
                pl.col("ensembl_gene_id").is_in(list_of_genes_in_region)
            )
            if annotations_genes.shape[0] == 0:
                continue
            # do a cross join between the subsetted annotations and the locus_id information
            cross_df = annotations_genes.select(
                [
                    pl.col("chr"),
                    pl.col("start"),
                    pl.col("end"),
                    pl.col("ensembl_gene_id"),
                    pl.col(score_col_dataset),
                ]
            ).join(
                locus_var_df.select(
                    [
                        pl.col("locus_id"),
                        pl.col("chr_hg19"),
                        pl.col("pos_hg19"),
                        pl.col(prob_col),
                    ]
                ),
                how="cross",
            )
            # filter for exact matches in region
            matches_df = cross_df.filter(
                (pl.col("chr") == pl.col("chr_hg19"))
                & (pl.col("pos_hg19") >= pl.col("start"))
                & (pl.col("pos_hg19") <= pl.col("end"))
            )
            if matches_df.shape[0] == 0:
                continue
            # calculate the scaled measurement (with probability or r2)
            # matches_df = matches_df.with_columns(
            #    (pl.col(score_col_dataset) * pl.col(prob_col)).alias(
            #        "scaled_measurement"
            #    )
            # )
            score_vals = matches_df[score_col_dataset].to_list()
            prob_vals = matches_df[prob_col].to_list()
            scaled_vals = [s * p for s, p in zip(score_vals, prob_vals)]

            matches_df = matches_df.with_columns(
                pl.Series("scaled_measurement", scaled_vals)
            )

            # keep keys and the measurement
            matches_df = matches_df.select(
                [
                    pl.col("locus_id"),
                    pl.col("ensembl_gene_id"),
                    pl.col("scaled_measurement"),
                ]
            )

            loci_annotations.append(matches_df)

        if loci_annotations:
            # concatenate the annotations from the different loci
            loci_annotations = pl.concat(loci_annotations)
        else:
            # typed empty frame
            loci_annotations = pl.DataFrame(
                {
                    "locus_id": pl.Series([], dtype=pl.String),
                    "ensembl_gene_id": pl.Series([], dtype=pl.String),
                    "scaled_measurement": pl.Series([], dtype=pl.Float64),
                }
            )
        return loci_annotations

    for method in ["CAGE", "PRO-seq_GRO-seq"]:
        dataset_name = f"HACER_{method}"
        sub_d_annotations = []
        for d in ["F5", "closest", "4D"]:
            # load HACER data from the per gene files
            dataset = _read_gene_parquet_files(
                os.path.join(config.paths.gene_regulation.hacer, method, d),
                genes_in_region["ensembl_gene_id"],
            )
            # if expected gene column missing, nothing to annotate for this subset
            if gene_col_dataset not in dataset.columns:
                sub_d_annotations.append(
                    pl.DataFrame(
                        {
                            "locus_id": pl.Series([], dtype=pl.String),
                            "ensembl_gene_id": pl.Series([], dtype=pl.String),
                            "scaled_measurement": pl.Series([], dtype=pl.Float64),
                        }
                    )
                )
                continue

            dataset = dataset.rename({gene_col_dataset: "ensembl_gene_id"})
            if "density" not in dataset.columns:
                dataset = dataset.with_columns(
                    pl.lit(1).alias(score_col_dataset).cast(pl.Float64)
                )
            dataset = dataset.with_columns(
                pl.col("start").cast(pl.Int64), pl.col("end").cast(pl.Int64)
            )
            logger.info(f" - - Loaded HACER {method} {d} data: {dataset.shape}")
            # get annotations per locus
            annotated_loci = get_subhacer_annotations(
                locus_groups_pl, dataset, genes_in_region
            )
            sub_d_annotations.append(annotated_loci)

        # concatenate the annotations from the different subdirectories
        if sub_d_annotations and any(df.shape[0] > 0 for df in sub_d_annotations):
            hacers = pl.concat(sub_d_annotations)
        else:
            # empty frame with the output columns
            hacers = pl.DataFrame(
                {
                    "locus_id": pl.Series([], dtype=pl.String),
                    "ensembl_gene_id": pl.Series([], dtype=pl.String),
                    "scaled_measurement": pl.Series([], dtype=pl.Float64),
                }
            )

        # aggregate the data with sum, max, mean
        if hacers.shape[0] == 0:
            results_df = pl.DataFrame(
                {
                    "locus_id": pl.Series([], dtype=pl.String),
                    "ensembl_gene_id": pl.Series([], dtype=pl.String),
                    f"{dataset_name}_sum": pl.Series([], dtype=pl.Float64),
                    f"{dataset_name}_max": pl.Series([], dtype=pl.Float64),
                    f"{dataset_name}_mean": pl.Series([], dtype=pl.Float64),
                }
            )
        else:
            results_df = hacers.group_by(["locus_id", "ensembl_gene_id"]).agg(
                pl.sum("scaled_measurement").alias(f"{dataset_name}_sum"),
                pl.max("scaled_measurement").alias(f"{dataset_name}_max"),
                pl.mean("scaled_measurement").alias(f"{dataset_name}_mean"),
            )

        # add the results to the hacer namespace
        hacer_results[dataset_name] = results_df

    return hacer_results["HACER_CAGE"], hacer_results["HACER_PRO-seq_GRO-seq"]


################################################
####           qtl colocalisation           ####
################################################


def get_GTEx_annotations(variant_data, genes_in_region, config, prob_col=None):
    """GTEx CAVIAR CLPP of the locus variants per candidate gene, sum and max."""

    gene_col_dataset = "GENE"  # ensembl_id
    # default to a constant weight of 1
    if not prob_col:
        prob_col = "prob_col"
        variant_data = variant_data.with_columns(
            pl.lit(1).alias(prob_col).cast(pl.Float64)
        )
    # load GTEx from the per gene files
    # the combined parquet adds `marker_name_hg38`, kept here for backwards compatibility
    # eqtl_gtex = pl.read_parquet(config.paths.eqtl.gtex_caviar_flames)
    eqtl_gtex = _read_gene_parquet_files(
        config.paths.eqtl.gtex_caviar_flames,
        genes_in_region["ensembl_gene_id"],
    )
    # the per gene files carry the raw `eQTL` field, the join key is derived from it
    eqtl_gtex = eqtl_gtex.with_columns(_marker_name_hg38(pl.col("eQTL")))
    eqtl_gtex = eqtl_gtex.select(
        [
            "GENE",
            "marker_name_hg38",
            "Probability",
        ]  # TODO include tissue weighting and TISSUE here
    )
    # split the variants by locus
    locus_groups = variant_data.partition_by("locus_id")

    # collect results per locus
    loci_annotations = []

    for locus_var_df in tqdm(locus_groups, desc="Processing loci for GTEx eQTLs"):
        locus_id = locus_var_df["locus_id"][0]

        # candidate genes of this locus
        gene_list = genes_in_region.filter(pl.col("locus_id") == locus_id)[
            "ensembl_gene_id"
        ]

        annotations_genes = eqtl_gtex.filter(pl.col(gene_col_dataset).is_in(gene_list))
        for gene_id in gene_list:
            # merge variant data with eQTLs
            coloc_variants_eqtl_gtex = locus_var_df.join(
                annotations_genes,
                left_on="marker_name_hg38",
                right_on="marker_name_hg38",
                how="inner",
            )
            coloc_variants_eqtl_gtex = coloc_variants_eqtl_gtex.unique()
            if coloc_variants_eqtl_gtex.shape[0] > 0:
                # calculate scaled probabilities
                coloc_variants_eqtl_gtex = coloc_variants_eqtl_gtex.with_columns(
                    (pl.col("Probability") * pl.col(prob_col)).alias("scaled_prob")
                )

                # aggregate CLPP scores by locus and gene
                results = pl.DataFrame(
                    {
                        "locus_id": [locus_id],
                        "ensembl_gene_id": [gene_id],
                        "CLPP_GTEX_eQTL_sum": [
                            coloc_variants_eqtl_gtex["scaled_prob"].sum()
                        ],
                        "CLPP_GTEX_eQTL_max": [
                            coloc_variants_eqtl_gtex["scaled_prob"].max()
                        ],
                    }
                )

                # add tissue-specific scores as a comment for future optionality
                # this would involve using tissue weights from magma_scores:
                # if config.use_tissue_weights:
                #     tissue_scores = merged_data.with_columns(
                #         pl.col("TISSUE").str.slice(0, 27).alias("tissue_short")
                #     )
                #     tissue_scores = tissue_scores.with_columns(
                #         pl.col("tissue_short").map_dict(magma_scores).alias("tissue_weight")
                #     )
                #     tissue_scores = tissue_scores.with_columns(
                #         (pl.col("Probability") * pl.col("tissue_weight") * pl.col(prob_col)).alias("tissue_scaled_prob")
                #     )
                #     results = results.with_columns([
                #         pl.lit(tissue_scores["tissue_scaled_prob"].sum()).alias("CLPP_GTEX_tissue_weighted_eQTL_sum"),
                #         pl.lit(tissue_scores["tissue_scaled_prob"].max()).alias("CLPP_GTEX_tissue_weighted_eQTL_max")
                #     ])

                loci_annotations.append(results)

    # stack the loci
    if loci_annotations:
        return pl.concat(loci_annotations)
    else:
        # empty frame with the output columns
        return pl.DataFrame(
            {
                "locus_id": [],
                "ensembl_gene_id": [],
                "CLPP_GTEX_eQTL_sum": [],
                "CLPP_GTEX_eQTL_max": [],
            }
        )


def get_eQTLgen_annotations(variant_data, genes_in_region, config, prob_col=None):
    """eQTLGen CLPP of the locus variants per candidate gene, sum and max."""
    gene_col_dataset = "ensembl_gene_id"  # ensembl_id

    # default to a constant weight of 1
    if not prob_col:
        prob_col = "prob_col"
        variant_data = variant_data.with_columns(
            pl.lit(1).alias(prob_col).cast(pl.Float64)
        )

    # load eQTLgen data from the per gene files
    # the combined parquet adds `marker_name_hg38` and the gene column, kept for compatibility
    # eqtl_gen = pl.read_parquet(config.paths.eqtl.eqtlgen)
    eqtl_gen = _read_gene_parquet_files(
        config.paths.eqtl.eqtlgen,
        genes_in_region["ensembl_gene_id"],
        gene_col=gene_col_dataset,
    )
    eqtl_gen = eqtl_gen.with_columns(_marker_name_hg38(pl.col("hg38")))
    eqtl_gen = eqtl_gen.select(
        [gene_col_dataset, "marker_name_hg38", "Causal_Post._Prob."]
    )

    # split the variants by locus
    locus_groups = variant_data.partition_by("locus_id")

    # collect results per locus
    loci_annotations = []

    for locus_var_df in tqdm(locus_groups, desc="Processing loci for eQTLgen"):
        # candidate genes of this locus
        gene_list = genes_in_region.filter(pl.col("locus_id") == locus_var_df["locus_id"][0])[
            "ensembl_gene_id"
        ]

        annotations_genes = eqtl_gen.filter(pl.col(gene_col_dataset).is_in(gene_list))

        for gene_id in gene_list:
            # merge variant data with eQTLs
            coloc_variants_eqtl_gen = locus_var_df.join(
                annotations_genes,
                left_on="marker_name_hg38",
                right_on="marker_name_hg38",
                how="inner",
            )
            coloc_variants_eqtl_gen = coloc_variants_eqtl_gen.unique()

            if coloc_variants_eqtl_gen.shape[0] > 0:
                # calculate scaled probabilities
                coloc_variants_eqtl_gen = coloc_variants_eqtl_gen.with_columns(
                    (pl.col("Causal_Post._Prob.") * pl.col(prob_col)).alias(
                        "scaled_prob"
                    )
                )

                # aggregate CLPP scores by locus and gene
                results = pl.DataFrame(
                    {
                        "locus_id": [locus_var_df["locus_id"][0]],
                        "ensembl_gene_id": [gene_id],
                        "CLPP_eQTLgen_sum": [
                            coloc_variants_eqtl_gen["scaled_prob"].sum()
                        ],
                        "CLPP_eQTLgen_max": [
                            coloc_variants_eqtl_gen["scaled_prob"].max()
                        ],
                    }
                )

                loci_annotations.append(results)

    # stack the loci
    if loci_annotations:
        return pl.concat(loci_annotations)
    else:
        # empty frame with the output columns
        return pl.DataFrame(
            {
                "locus_id": [],
                "ensembl_gene_id": [],
                "CLPP_eQTLgen_sum": [],
                "CLPP_eQTLgen_max": [],
            }
        )


def get_QTL_catalog_annotations(variant_data, genes_in_region, config, prob_col=None):
    """eQTL Catalogue eQTL and rQTL CLPP of the locus variants per candidate gene."""
    gene_col_dataset = "ensg"  # ensembl_id

    # default to a constant weight of 1
    if not prob_col:
        prob_col = "prob_col"
        variant_data = variant_data.with_columns(
            pl.lit(1).alias(prob_col).cast(pl.Float64)
        )

    results = {}

    # process both eQTL and rQTL data
    for qtl_type in ["eQTL", "rQTL"]:
        # load corresponding QTL catalog data from the per gene files
        # the combined parquets add `marker_name_hg38`, kept here for backwards compatibility
        # if qtl_type == "eQTL":
        #     qtl_catalog = pl.read_parquet(config.paths.eqtl.eqtl_catalogue_eqtl)
        # else:
        #     qtl_catalog = pl.read_parquet(config.paths.eqtl.eqtl_catalogue_rqtl)
        if qtl_type == "eQTL":
            qtl_directory = config.paths.eqtl.eqtl_catalogue_eqtl
        else:
            qtl_directory = config.paths.eqtl.eqtl_catalogue_rqtl

        qtl_catalog = _read_gene_parquet_files(
            qtl_directory, genes_in_region["ensembl_gene_id"]
        )
        # the per gene files carry the raw `QTL_b38` field, the join key is derived from it
        qtl_catalog = qtl_catalog.with_columns(_marker_name_hg38(pl.col("QTL_b38")))
        qtl_catalog = qtl_catalog.select(
            [gene_col_dataset, "marker_name_hg38", "max", "mean"]
        )

        # split the variants by locus
        locus_groups = variant_data.partition_by("locus_id")

        loci_annotations = []

        for locus_var_df in tqdm(
            locus_groups, desc=f"Processing loci for QTL catalog {qtl_type}"
        ):

            # candidate genes of this locus
            gene_list = genes_in_region.filter(pl.col("locus_id") == locus_var_df["locus_id"][0])[
                "ensembl_gene_id"
            ]

            annotations_genes = qtl_catalog.filter(
                pl.col(gene_col_dataset).is_in(gene_list)
            )

            for gene_id in gene_list:
                # merge variant data with QTLs
                coloc_variants_qtl = locus_var_df.join(
                    annotations_genes,
                    left_on="marker_name_hg38",
                    right_on="marker_name_hg38",
                    how="inner",
                )
                coloc_variants_qtl = coloc_variants_qtl.unique()

                if coloc_variants_qtl.shape[0] > 0:
                    # calculate scaled probabilities using 'max' and 'mean' columns
                    coloc_variants_qtl = coloc_variants_qtl.with_columns(
                        [
                            (pl.col("max") * pl.col(prob_col)).alias("max_scaled"),
                            (pl.col("mean") * pl.col(prob_col)).alias("mean_scaled"),
                        ]
                    )

                    # aggregate scores
                    results_df = pl.DataFrame(
                        {
                            "locus_id": [locus_var_df["locus_id"][0]],
                            "ensembl_gene_id": [gene_id],
                            f"max_CLPP_{qtl_type}_eQTLCatalog_sum": [
                                coloc_variants_qtl["max_scaled"].sum()
                            ],
                            f"max_CLPP_{qtl_type}_eQTLCatalog_max": [
                                coloc_variants_qtl["max_scaled"].max()
                            ],
                            f"mean_CLPP_{qtl_type}_eQTLCatalog_sum": [
                                coloc_variants_qtl["mean_scaled"].sum()
                            ],
                            f"mean_CLPP_{qtl_type}_eQTLCatalog_max": [
                                coloc_variants_qtl["mean_scaled"].max()
                            ],
                        }
                    )

                    loci_annotations.append(results_df)

        # stack the loci for this QTL type
        if loci_annotations:
            results[qtl_type] = pl.concat(loci_annotations)
        else:
            # empty frame with the output columns
            results[qtl_type] = pl.DataFrame(
                {
                    "locus_id": [],
                    "ensembl_gene_id": [],
                    f"max_CLPP_{qtl_type}_eQTLCatalog_sum": [],
                    f"max_CLPP_{qtl_type}_eQTLCatalog_max": [],
                    f"mean_CLPP_{qtl_type}_eQTLCatalog_sum": [],
                    f"mean_CLPP_{qtl_type}_eQTLCatalog_max": [],
                }
            )

    return results["eQTL"], results["rQTL"]


################################################
####           chromatin contacts           ####
################################################


def get_hic_annotations(variant_data, genes_in_region, config, prob_col=None):
    """Jung and Javierre Hi-C contact scores per (locus, gene), returned as (jung, javierre)."""

    def process_hic_dataset(locus_df, hic_data, gene_list, measurement_col, prob_col):
        """Weighted Hi-C measurements of the locus variants inside each gene's contact intervals, as dicts."""
        # filter HiC data for the genes of interest
        annotations_genes = hic_data.filter(pl.col("ENSG").is_in(gene_list))

        if annotations_genes.shape[0] == 0:
            return []

        # create a cross join between HiC data and locus information
        cross_df = annotations_genes.select(
            [
                pl.col("chr"),
                pl.col("start"),
                pl.col("end"),
                pl.col("ENSG").alias("ensembl_gene_id"),
                pl.col(measurement_col),
            ]
        ).join(
            locus_df.select(
                [
                    pl.col("locus_id"),
                    pl.col("chr_hg19"),
                    pl.col("pos_hg19"),
                    pl.col(prob_col),
                ]
            ),
            how="cross",
        )

        # filter for exact matches in region
        matches_df = cross_df.filter(
            (pl.col("chr") == pl.col("chr_hg19"))
            & (pl.col("pos_hg19") >= pl.col("start"))
            & (pl.col("pos_hg19") <= pl.col("end"))
        )

        if matches_df.shape[0] == 0:
            return []

        # calculate the scaled measurement
        matches_df = matches_df.with_columns(
            (pl.col(measurement_col) * pl.col(prob_col)).alias("scaled_measurement")
        )

        # convert to list of dictionaries
        result = matches_df.select(
            [
                pl.col("locus_id"),
                pl.col("ensembl_gene_id"),
                pl.col("scaled_measurement"),
            ]
        ).to_dicts()

        return result

    # default to a constant weight of 1
    if not prob_col:
        prob_col = "prob_col"
        variant_data = variant_data.with_columns(
            pl.lit(1).alias(prob_col).cast(pl.Float64)
        )

    # load Javierre HiC data from the per gene files
    # the combined parquet adds a hg38 liftover, kept here for backwards compatibility
    # javierre_data = pl.read_parquet(
    #     config.paths.gene_regulation.hic_annotation.javierre_hic_data
    # )
    javierre_data = _read_gene_parquet_files(
        config.paths.gene_regulation.hic_annotation.javierre_hic_data,
        genes_in_region["ensembl_gene_id"],
    )
    # the hg38 columns were lifted from the other end, so hg19 matching uses `oe` directly
    javierre_data = javierre_data.with_columns(
        pl.concat_str([pl.lit("chr"), pl.col("oeChr").cast(pl.String)]).alias("chr"),
        pl.col("oeStart").alias("start"),
        pl.col("oeEnd").alias("end"),
        pl.sum_horizontal(JAVIERRE_CELL_TYPES).alias("sum"),
        pl.max_horizontal(JAVIERRE_CELL_TYPES).alias("max"),
    )

    # load Jung HiC data from the per gene files
    # jung_data = pl.read_parquet(
    #     config.paths.gene_regulation.hic_annotation.jung_hic_data
    # )
    jung_data = _read_gene_parquet_files(
        config.paths.gene_regulation.hic_annotation.jung_hic_data,
        genes_in_region["ensembl_gene_id"],
    )
    # `target` holds the hg19 interval as chr.start.end, and `p` is a copy of `dist_pvalue`
    jung_data = jung_data.with_columns(
        pl.col("target").str.split(".").list.get(1).cast(pl.Int64).alias("start"),
        pl.col("target").str.split(".").list.get(2).cast(pl.Int64).alias("end"),
        pl.col("dist_pvalue").alias("p"),
    )

    # split the variants by locus
    locus_groups_pl = variant_data.partition_by("locus_id")

    # create empty lists to store annotation results
    jung_annotations = []
    javierre_sum_annotations = []
    javierre_max_annotations = []

    # iterate through each locus
    for locus_var_df in tqdm(locus_groups_pl, desc="Processing loci for HiC annotation"):
        # candidate genes of this locus
        list_of_genes_in_region = genes_in_region.filter(
            pl.col("locus_id") == locus_var_df["locus_id"][0]
        )["ensembl_gene_id"]

        # process Javierre data
        javierre_sum_annotations.extend(
            process_hic_dataset(
                locus_var_df, javierre_data, list_of_genes_in_region, "sum", prob_col
            )
        )

        javierre_max_annotations.extend(
            process_hic_dataset(
                locus_var_df, javierre_data, list_of_genes_in_region, "max", prob_col
            )
        )

        # process Jung data
        jung_annotations.extend(
            process_hic_dataset(
                locus_var_df, jung_data, list_of_genes_in_region, "p", prob_col
            )
        )

    # create dataframes from the collected annotations
    if jung_annotations:
        jung_results = pl.DataFrame(jung_annotations)
        jung_results = jung_results.group_by(["locus_id", "ensembl_gene_id"]).agg(
            pl.sum("scaled_measurement").alias("Jung_HiC_sum"),
            pl.max("scaled_measurement").alias("Jung_HiC_max"),
            pl.mean("scaled_measurement").alias("Jung_HiC_mean"),
        )
    else:
        jung_results = pl.DataFrame(
            {
                "locus_id": [],
                "ensembl_gene_id": [],
                "Jung_HiC_sum": [],
                "Jung_HiC_max": [],
                "Jung_HiC_mean": [],
            }
        )

    # combine Javierre sum and max results
    if javierre_sum_annotations and javierre_max_annotations:
        javierre_sum_results = pl.DataFrame(javierre_sum_annotations)
        javierre_sum_results = javierre_sum_results.group_by(
            ["locus_id", "ensembl_gene_id"]
        ).agg(
            pl.sum("scaled_measurement").alias("Javierre_HiC_sum"),
            pl.mean("scaled_measurement").alias("Javierre_HiC_mean"),
        )

        javierre_max_results = pl.DataFrame(javierre_max_annotations)
        javierre_max_results = javierre_max_results.group_by(
            ["locus_id", "ensembl_gene_id"]
        ).agg(
            pl.max("scaled_measurement").alias("Javierre_HiC_max"),
        )

        javierre_results = javierre_sum_results.join(
            javierre_max_results, on=["locus_id", "ensembl_gene_id"], how="outer"
        )
    else:
        javierre_results = pl.DataFrame(
            {
                "locus_id": [],
                "ensembl_gene_id": [],
                "Javierre_HiC_sum": [],
                "Javierre_HiC_max": [],
                "Javierre_HiC_mean": [],
            }
        )

    return jung_results, javierre_results


################################################
####        promoters and enhancers         ####
################################################


def get_Promoters_annotations(variant_data, genes_in_region, config, prob_col=None):
    """Summed and max weight of the locus variants inside each candidate gene's promoters."""
    # default to a constant weight of 1
    if not prob_col:
        prob_col = "prob_col"
        variant_data = variant_data.with_columns(
            pl.lit(1).alias(prob_col).cast(pl.Float64)
        )

    # load promoter data
    promoters = pl.read_parquet(config.paths.gene_regulation.promoter)
    # the FLAMES file ships bare chromosome names, `chr_hg38` carries the `chr` prefix
    promoters = promoters.with_columns(
        pl.concat_str(
            [pl.lit("chr"), pl.col("chromosome").str.replace("^chr", "")]
        ).alias("chromosome")
    )

    # split the variants by locus
    locus_groups = variant_data.partition_by("locus_id")

    # collect results per locus
    loci_annotations = []

    for locus_var_df in tqdm(locus_groups, desc="Processing loci for promoter annotations"):
        locus_id = locus_var_df["locus_id"][0]

        # candidate genes of this locus
        gene_list = genes_in_region.filter(pl.col("locus_id") == locus_id)[
            "ensembl_gene_id"
        ]

        # position columns based on genome build
        chr_col = "chr_hg38"
        pos_col = "pos_hg38"

        for gene_id in gene_list:
            # filter promoters for this gene
            gene_promoters = promoters.filter(pl.col("ensg") == gene_id)

            if gene_promoters.shape[0] == 0:
                # no promoters found for this gene
                continue

            # cross join locus variants with gene promoters
            cross_df = gene_promoters.select(
                [pl.col("chromosome"), pl.col("start"), pl.col("end"), pl.col("ensg")]
            ).join(
                locus_var_df.select(
                    [
                        pl.col("locus_id"),
                        pl.col(chr_col),
                        pl.col(pos_col),
                        pl.col(prob_col),
                    ]
                ),
                how="cross",
            )

            # filter for variants in promoter regions
            matches_df = cross_df.filter(
                (pl.col("chr_hg38") == pl.col("chromosome"))
                & (pl.col(pos_col) >= pl.col("start"))
                & (pl.col(pos_col) <= pl.col("end"))
            )

            if matches_df.shape[0] > 0:
                # aggregate results for this gene
                results_df = pl.DataFrame(
                    {
                        "locus_id": [locus_id],
                        "ensembl_gene_id": [gene_id],
                        "Promoter_sum": [matches_df[prob_col].sum()],
                        "Promoter_max": [matches_df[prob_col].max()],
                    }
                )

                loci_annotations.append(results_df)

    # stack the loci
    if loci_annotations:
        return pl.concat(loci_annotations)
    else:
        # empty frame with the output columns
        return pl.DataFrame(
            {
                "locus_id": [],
                "ensembl_gene_id": [],
                "Promoter_sum": [],
                "Promoter_max": [],
            }
        )


def get_F5_annotations(variant_data, genes_in_region, config, prob_col=None):
    """FANTOM5 enhancer-gene scores of the locus variants per candidate gene."""
    dataset_name = "F5_ep"

    # default to a constant weight of 1
    if not prob_col:
        prob_col = "prob_col"
        variant_data = variant_data.with_columns(
            pl.lit(1).alias(prob_col).cast(pl.Float64)
        )

    # load FANTOM5 data from the per gene files
    # the combined parquet is a self-made derivative, kept here for backwards compatibility
    # f5_data = pl.read_parquet(config.paths.gene_regulation.fantom5)
    f5_data = _read_gene_parquet_files(
        config.paths.gene_regulation.fantom5,
        genes_in_region["ensembl_gene_id"],
    )
    f5_data = f5_data.rename({"ENSG_promoter": "ensembl_gene_id"})

    logger.info(f"Loaded FANTOM5 data: {f5_data.shape}")

    # split the variants by locus
    locus_groups_pl = variant_data.partition_by("locus_id")

    # initialize list for annotations
    loci_annotations = []

    for locus_var_df in tqdm(locus_groups_pl, desc="Processing loci for FANTOM5"):
        list_of_genes_in_region = genes_in_region.filter(
            pl.col("locus_id") == locus_var_df["locus_id"][0]
        )["ensembl_gene_id"]

        annotations_genes = f5_data.filter(
            pl.col("ensembl_gene_id").is_in(list_of_genes_in_region)
        )
        # rename columns for consistency
        if annotations_genes.shape[0] == 0:
            continue

        # extract chr, start, end from enhancer field if needed
        if "enhancer" in annotations_genes.columns and not all(
            col in annotations_genes.columns for col in ["chr", "start", "end"]
        ):
            annotations_genes = annotations_genes.with_columns(
                [
                    pl.col("enhancer").str.extract(r"^([^:]+):", 1).alias("chr"),
                    pl.col("enhancer")
                    .str.extract(r":(\d+)-", 1)
                    .cast(pl.Int64)
                    .alias("start"),
                    pl.col("enhancer")
                    .str.extract(r"-(\d+)$", 1)
                    .cast(pl.Int64)
                    .alias("end"),
                ]
            )

        # do a cross join between the annotations and locus variants
        # we check if the variants are in enhancer regions for a given gene
        cross_df = annotations_genes.select(
            [
                pl.col("chr"),
                pl.col("start"),
                pl.col("end"),
                pl.col("ensembl_gene_id"),
                pl.col("correlation"),
            ]
        ).join(
            locus_var_df.select(
                [
                    pl.col("locus_id"),
                    pl.col("chr_hg19"),
                    pl.col("pos_hg19"),
                    pl.col(prob_col),
                ]
            ),
            how="cross",
        )

        # filter for exact matches in region
        matches_df = cross_df.filter(
            (pl.col("chr") == pl.col("chr_hg19"))
            & (pl.col("pos_hg19") >= pl.col("start"))
            & (pl.col("pos_hg19") <= pl.col("end"))
        )

        if matches_df.shape[0] == 0:
            continue

        # calculate scaled correlation
        matches_df = matches_df.with_columns(
            (pl.col("correlation") * pl.col(prob_col)).alias("scaled_correlation")
        )

        # aggregate results
        results_df = matches_df.group_by(["locus_id", "ensembl_gene_id"]).agg(
            pl.sum("scaled_correlation").alias(f"{dataset_name}_sum"),
            pl.max("scaled_correlation").alias(f"{dataset_name}_max"),
        )

        loci_annotations.append(results_df)

    # combine results
    if loci_annotations:
        return pl.concat(loci_annotations)
    else:
        return pl.DataFrame(
            {
                "locus_id": [],
                "ensembl_gene_id": [],
                f"{dataset_name}_sum": [],
                f"{dataset_name}_max": [],
            }
        )


def get_GeneHancer_annotations(variant_data, genes_in_region, config, prob_col=None):
    """GeneHancer enhancer-gene scores of the locus variants per candidate gene."""
    dataset_name = "GeneHancer"

    # default to a constant weight of 1
    if not prob_col:
        prob_col = "prob_col"
        variant_data = variant_data.with_columns(
            pl.lit(1).alias(prob_col).cast(pl.Float64)
        )

    # load GeneHancer data from the per gene files
    # the combined parquet is a self-made derivative, kept here for backwards compatibility
    # genehancer_data = pl.read_parquet(config.paths.gene_regulation.genehancer)
    genehancer_data = _read_gene_parquet_files(
        config.paths.gene_regulation.genehancer,
        genes_in_region["ensembl_gene_id"],
    )
    genehancer_data = genehancer_data.rename({"gene": "ensembl_gene_id"})
    logger.info(f"Loaded GeneHancer data: {genehancer_data.shape}")

    # split the variants by locus
    locus_groups_pl = variant_data.partition_by("locus_id")

    # initialize list for annotations
    loci_annotations = []

    for locus_id in tqdm(locus_groups_pl, desc="Processing loci for GeneHancer"):
        list_of_genes_in_region = genes_in_region.filter(
            pl.col("locus_id") == locus_id["locus_id"][0]
        )["ensembl_gene_id"]

        annotations_genes = genehancer_data.filter(
            pl.col("ensembl_gene_id").is_in(list_of_genes_in_region)
        )

        if annotations_genes.shape[0] == 0:
            continue

        # do a cross join between the annotations and locus variants
        cross_df = annotations_genes.select(
            [
                pl.col("chr"),
                pl.col("start"),
                pl.col("end"),
                pl.col("ensembl_gene_id"),
                pl.col("score"),
            ]
        ).join(
            locus_id.select(
                [
                    pl.col("locus_id"),
                    pl.col("chr_hg19"),
                    pl.col("pos_hg19"),
                    pl.col(prob_col),
                ]
            ),
            how="cross",
        )

        # filter for exact matches in region
        matches_df = cross_df.filter(
            (pl.col("chr") == pl.col("chr_hg19"))
            & (pl.col("pos_hg19") >= pl.col("start"))
            & (pl.col("pos_hg19") <= pl.col("end"))
        )

        if matches_df.shape[0] == 0:
            continue

        # calculate scaled score
        matches_df = matches_df.with_columns(
            (pl.col("score") * pl.col(prob_col)).alias("scaled_score")
        )

        # aggregate results
        results_df = matches_df.group_by(["locus_id", "ensembl_gene_id"]).agg(
            pl.sum("scaled_score").alias(f"{dataset_name}_sum"),
            pl.max("scaled_score").alias(f"{dataset_name}_max"),
        )

        loci_annotations.append(results_df)

    # combine results
    if loci_annotations:
        return pl.concat(loci_annotations)
    else:
        return pl.DataFrame(
            {
                "locus_id": [],
                "ensembl_gene_id": [],
                f"{dataset_name}_sum": [],
                f"{dataset_name}_max": [],
            }
        )


def get_RoadmapEpi_annotations(variant_data, genes_in_region, config, prob_col=None):
    """Roadmap Epigenomics enhancer-gene scores of the locus variants per candidate gene."""
    dataset_name = "RoadmapEpi"

    # default to a constant weight of 1
    if not prob_col:
        prob_col = "prob_col"
        variant_data = variant_data.with_columns(
            pl.lit(1).alias(prob_col).cast(pl.Float64)
        )

    # load Roadmap data from the per gene files
    # the combined parquet is a self-made derivative, kept here for backwards compatibility
    # roadmap_data = pl.read_parquet(config.paths.gene_regulation.roadmap)
    roadmap_data = _read_gene_parquet_files(
        config.paths.gene_regulation.roadmap,
        genes_in_region["ensembl_gene_id"],
    )
    roadmap_data = roadmap_data.rename({"ENSG": "ensembl_gene_id"})
    logger.info(f"Loaded Roadmap Epigenomics data: {roadmap_data.shape}")

    # split the variants by locus
    locus_groups_pl = variant_data.partition_by("locus_id")

    # initialize list for annotations
    loci_annotations = []

    for locus_var_df in tqdm(
        locus_groups_pl, desc="Processing loci for Roadmap Epigenomics"
    ):
        list_of_genes_in_region = genes_in_region.filter(
            pl.col("locus_id") == locus_var_df["locus_id"][0]
        )["ensembl_gene_id"]

        annotations_genes = roadmap_data.filter(
            pl.col("ensembl_gene_id").is_in(list_of_genes_in_region)
        )

        if annotations_genes.shape[0] == 0:
            continue

        # do a cross join between the annotations and locus variants
        cross_df = annotations_genes.select(
            [
                pl.col("chr"),
                pl.col("start"),
                pl.col("end"),
                pl.col("ensembl_gene_id"),
                pl.col("max"),
            ]
        ).join(
            locus_var_df.select(
                [
                    pl.col("locus_id"),
                    pl.col("chr_hg19"),
                    pl.col("pos_hg19"),
                    pl.col(prob_col),
                ]
            ),
            how="cross",
        )

        # filter for exact matches in region
        matches_df = cross_df.filter(
            (pl.col("chr") == pl.col("chr_hg19"))
            & (pl.col("pos_hg19") >= pl.col("start"))
            & (pl.col("pos_hg19") <= pl.col("end"))
        )

        if matches_df.shape[0] == 0:
            continue

        # calculate scaled max score
        matches_df = matches_df.with_columns(
            (pl.col("max") * pl.col(prob_col)).alias("scaled_max")
        )

        # aggregate results
        results_df = matches_df.group_by(["locus_id", "ensembl_gene_id"]).agg(
            pl.sum("scaled_max").alias(f"{dataset_name}_sum"),
            pl.max("scaled_max").alias(f"{dataset_name}_max"),
        )

        loci_annotations.append(results_df)

    # combine results
    if loci_annotations:
        return pl.concat(loci_annotations)
    else:
        return pl.DataFrame(
            {
                "locus_id": [],
                "ensembl_gene_id": [],
                f"{dataset_name}_sum": [],
                f"{dataset_name}_max": [],
            }
        )


def get_EpiMap_annotations(variant_data, genes_in_region, config, prob_col=None):
    """EpiMap enhancer-gene scores of the locus variants per candidate gene."""
    dataset_name = "EpiMap"

    # default to a constant weight of 1
    if not prob_col:
        prob_col = "prob_col"
        variant_data = variant_data.with_columns(
            pl.lit(1).alias(prob_col).cast(pl.Float64)
        )

    # load EpiMap data from the per gene files
    # the combined parquet is a self-made derivative, kept here for backwards compatibility
    # epimap_data = pl.read_parquet(config.paths.gene_regulation.epimap)
    epimap_data = _read_gene_parquet_files(
        config.paths.gene_regulation.epimap,
        genes_in_region["ensembl_gene_id"],
    )
    epimap_data = epimap_data.rename({"ensg": "ensembl_gene_id"})
    logger.info(f"Loaded EpiMap data: {epimap_data.shape}")

    # split the variants by locus
    locus_groups_pl = variant_data.partition_by("locus_id")

    # initialize list for annotations
    loci_annotations = []

    for locus_id in tqdm(locus_groups_pl, desc="Processing loci for EpiMap"):
        list_of_genes_in_region = genes_in_region.filter(
            pl.col("locus_id") == locus_id["locus_id"][0]
        )["ensembl_gene_id"]

        annotations_genes = epimap_data.filter(
            pl.col("ensembl_gene_id").is_in(list_of_genes_in_region)
        )

        if annotations_genes.shape[0] == 0:
            continue

        # do a cross join between the annotations and locus variants
        cross_df = annotations_genes.select(
            [
                pl.col("chr"),
                pl.col("start"),
                pl.col("end"),
                pl.col("ensembl_gene_id"),
                pl.col("corr"),
            ]
        ).join(
            locus_id.select(
                [
                    pl.col("locus_id"),
                    pl.col("chr_hg19"),
                    pl.col("pos_hg19"),
                    pl.col(prob_col),
                ]
            ),
            how="cross",
        )

        # filter for exact matches in region
        matches_df = cross_df.filter(
            (pl.col("chr") == pl.col("chr_hg19"))
            & (pl.col("pos_hg19") >= pl.col("start"))
            & (pl.col("pos_hg19") <= pl.col("end"))
        )

        if matches_df.shape[0] == 0:
            continue

        # calculate scaled correlation
        matches_df = matches_df.with_columns(
            (pl.col("corr") * pl.col(prob_col)).alias("scaled_corr")
        )

        # aggregate results
        results_df = matches_df.group_by(["locus_id", "ensembl_gene_id"]).agg(
            pl.sum("scaled_corr").alias(f"{dataset_name}_sum"),
            pl.max("scaled_corr").alias(f"{dataset_name}_max"),
        )

        loci_annotations.append(results_df)

    # combine results
    if loci_annotations:
        return pl.concat(loci_annotations)
    else:
        return pl.DataFrame(
            {
                "locus_id": [],
                "ensembl_gene_id": [],
                f"{dataset_name}_sum": [],
                f"{dataset_name}_max": [],
            }
        )


def get_CICERO_annotations(variant_data, genes_in_region, config, prob_col=None):
    """CICERO co-accessibility scores of the locus variants per candidate gene."""
    dataset_name = "cicero"

    # default to a constant weight of 1
    if not prob_col:
        prob_col = "prob_col"
        variant_data = variant_data.with_columns(
            pl.lit(1).alias(prob_col).cast(pl.Float64)
        )

    # load CICERO data from the per gene files
    # the combined parquet is a self-made derivative, kept here for backwards compatibility
    # cicero_data = pl.read_parquet(config.paths.gene_regulation.cicero)
    cicero_data = _read_gene_parquet_files(
        config.paths.gene_regulation.cicero,
        genes_in_region["ensembl_gene_id"],
    )
    cicero_data = cicero_data.rename({"ensg": "ensembl_gene_id"})
    logger.info(f"Loaded CICERO data: {cicero_data.shape}")

    # ensure chromosome format is consistent

    # split the variants by locus
    locus_groups_pl = variant_data.partition_by("locus_id")

    # initialize list for annotations
    loci_annotations = []

    for locus_var_df in tqdm(locus_groups_pl, desc="Processing loci for CICERO"):
        list_of_genes_in_region = genes_in_region.filter(
            pl.col("locus_id") == locus_var_df["locus_id"][0]
        )["ensembl_gene_id"]

        annotations_genes = cicero_data.filter(
            pl.col("ensembl_gene_id").is_in(list_of_genes_in_region)
        )

        if annotations_genes.shape[0] == 0:
            continue

        # do a cross join between the annotations and locus variants
        cross_df = annotations_genes.select(
            [
                pl.col("chr"),
                pl.col("start"),
                pl.col("end"),
                pl.col("ensembl_gene_id"),
                pl.col("cor"),  # Correlation score in CICERO
            ]
        ).join(
            locus_var_df.select(
                [
                    pl.col("locus_id"),
                    pl.col(
                        "chr_hg19"
                    ),  # Using hg19 coordinates as CICERO is based on hg19
                    pl.col("pos_hg19"),
                    pl.col(prob_col),
                ]
            ),
            how="cross",
        )

        # filter for exact matches in region
        matches_df = cross_df.filter(
            (pl.col("chr") == pl.col("chr_hg19"))
            & (pl.col("pos_hg19") >= pl.col("start"))
            & (pl.col("pos_hg19") <= pl.col("end"))
        )

        if matches_df.shape[0] == 0:
            continue

        # calculate scaled correlation
        matches_df = matches_df.with_columns(
            (pl.col("cor") * pl.col(prob_col)).alias("scaled_cor")
        )

        # aggregate results
        results_df = matches_df.group_by(["locus_id", "ensembl_gene_id"]).agg(
            pl.sum("scaled_cor").alias(f"{dataset_name}_sum"),
            pl.max("scaled_cor").alias(f"{dataset_name}_max"),
            pl.mean("scaled_cor").alias(f"{dataset_name}_mean"),
        )

        loci_annotations.append(results_df)

    # combine results
    if loci_annotations:
        return pl.concat(loci_annotations)
    else:
        return pl.DataFrame(
            {
                "locus_id": [],
                "ensembl_gene_id": [],
                f"{dataset_name}_sum": [],
                f"{dataset_name}_max": [],
                f"{dataset_name}_mean": [],
            }
        )


################################################
####              entry point               ####
################################################


def get_flames_annotations(
    variant_data,
    genes_in_region,
    config,
    prob_col="scaling_factor",
):
    # weight each variant by prob_col, PIP for credible sets, r2 otherwise, 1.0 for leads
    # fall back to unweighted when the column is absent
    """Run every FLAMES annotation, each weighted by `prob_col`, and return them by name."""
    if prob_col and prob_col not in variant_data.columns:
        logger.warning(
            f"'{prob_col}' not found in variant_data; FLAMES annotations will be unweighted."
        )
        prob_col = None

    # get HACER annotations
    hacer_cage, hacer_pro_seq = get_HACER_annotations(
        variant_data, genes_in_region, config, prob_col=prob_col
    )

    # get ABC EP annotations
    abc_ep_annotation = get_ABC_EP_annotations(
        variant_data, genes_in_region, config, prob_col=prob_col
    )

    # get ABC CRISPR annotations
    abc_crispr_annotation = get_ABC_CRISPR_annotations(
        variant_data, genes_in_region, config, prob_col=prob_col
    )

    # get CICERO annotations
    cicero_annotations = get_CICERO_annotations(
        variant_data, genes_in_region, config, prob_col=prob_col
    )
    # get FANTOM5 annotations
    f5_annotations = get_F5_annotations(
        variant_data, genes_in_region, config, prob_col=prob_col
    )

    # get GeneHancer annotations
    genehancer_annotations = get_GeneHancer_annotations(
        variant_data, genes_in_region, config, prob_col=prob_col
    )

    # get Roadmap Epigenomics annotations
    roadmap_annotations = get_RoadmapEpi_annotations(
        variant_data, genes_in_region, config, prob_col=prob_col
    )

    # get EpiMap annotations
    epimap_annotations = get_EpiMap_annotations(
        variant_data, genes_in_region, config, prob_col=prob_col
    )

    # get promoter region annotations
    promoter_annotations = get_Promoters_annotations(
        variant_data, genes_in_region, config, prob_col=prob_col
    )

    # get eQTL GTEx CAVIAR finemapping annotations
    gtex_caviar_flames = get_GTEx_annotations(
        variant_data, genes_in_region, config, prob_col=prob_col
    )

    # get eQTLgen annotations
    eQTLgen_annotations = get_eQTLgen_annotations(
        variant_data, genes_in_region, config, prob_col=prob_col
    )

    # get QTL catalog annotations
    qtl_catalog_eQTL, qtl_catalog_rQTL = get_QTL_catalog_annotations(
        variant_data, genes_in_region, config, prob_col=prob_col
    )

    # get HIC annotations
    jung_hic, javierre_hic = get_hic_annotations(
        variant_data, genes_in_region, config, prob_col=prob_col
    )

    return {
        "promotor_annotation": promoter_annotations,
        "abc_ep_annotation": abc_ep_annotation,
        "abc_crispr_annotation": abc_crispr_annotation,
        "jung_hic": jung_hic,
        "javierre_hic": javierre_hic,
        "hacer_cage": hacer_cage,
        "hacer_pro_seq": hacer_pro_seq,
        "gtex_caviar_flames": gtex_caviar_flames,
        "eqtlgen_annotations": eQTLgen_annotations,
        "qtl_catalog_eQTL": qtl_catalog_eQTL,
        "qtl_catalog_rQTL": qtl_catalog_rQTL,
        "f5_annotations": f5_annotations,
        "genehancer_annotations": genehancer_annotations,
        "roadmap_annotations": roadmap_annotations,
        "epimap_annotations": epimap_annotations,
        "cicero_annotations": cicero_annotations,
    }
