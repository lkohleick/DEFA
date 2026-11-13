################################################
#### Candidate genes: window and distances  ####
################################################
"""Candidate genes around each lead and its weighted centroid, with lead and centroid distances."""

import logging
from pathlib import Path
import tqdm
import polars as pl
from pyensembl import EnsemblRelease

# module logger and repo root
logger = logging.getLogger(__name__)
project_root = Path(__file__).parent.parent.parent


################################################
####             locus regions              ####
################################################


def compute_regions(variant_data: pl.DataFrame) -> pl.DataFrame:
    """Per locus: min start, max end, variant count and the `scaling_factor`-weighted centroid."""

    # return the empty schema when there is no data
    if variant_data.is_empty():
        return pl.DataFrame(
            schema={
                "locus_id": pl.String,
                "region_start": pl.Int64,
                "region_end": pl.Int64,
                "region_centroid_weighted": pl.Int64,
                "num_variants": pl.UInt32,
            }
        )

    # cast positions and weights to numeric
    df = variant_data.with_columns(
        [pl.col("start_hg38").cast(pl.Int64), pl.col("scaling_factor").cast(pl.Float64)]
    )

    # --> aggregate per locus <-- #
    return df.group_by("locus_id").agg(
        [
            pl.col("start_hg38").min().alias("region_start"),
            pl.col("end_hg38").max().alias("region_end"),
            pl.count().alias("num_variants"),
            # weighted centroid, sum(position * weight) / sum(weight)
            (
                (pl.col("start_hg38") * pl.col("scaling_factor")).sum()
                / pl.col("scaling_factor").sum()
            )
            .cast(pl.Int64)
            .alias("region_centroid_weighted"),
        ]
    )


################################################
####            candidate genes             ####
################################################


def get_genes_in_region(
    variant_data: pl.DataFrame,
    config,
) -> pl.DataFrame:
    """Genes near the lead and the centroid, with footprint and TSS distances from both."""
    logger.info("Finding genes in region ...")

    # set up Ensembl and the window
    ensembl = EnsemblRelease(config.parameters.ensembl_release)
    biotype = config.parameters.biotypes_of_interest
    window_size = config.parameters.interval_kb * 1000

    # --> weighted centroids <-- #
    logger.info("Computing weighted centroids...")
    region_centroids = compute_regions(variant_data)

    # --> leads <-- #
    # the lead is the variant whose variant_id equals its locus_id
    lead_variants = variant_data.filter(
        pl.col("locus_id") == pl.col("variant_id")
    ).select(
        [
            "locus_id",
            "chr_hg38_ensembl",
            pl.col("start_hg38").cast(pl.Int64).alias("lead_pos"),
            "marker_name_hg38",
        ]
    )

    # join the regions back onto the leads
    locus_info_df = lead_variants.join(region_centroids, on="locus_id", how="left")

    result_rows = []

    # --> query genes per locus <-- #
    for row in tqdm.tqdm(locus_info_df.to_dicts(),desc="Getting genes for regions"):
        locus_id = row["locus_id"]
        chrom = row["chr_hg38_ensembl"]
        lead_pos = row["lead_pos"]
        marker_name = row["marker_name_hg38"]

        # region info, the lead position stands in when the centroid is missing
        centroid = row.get("region_centroid_weighted") or lead_pos
        r_start = row.get("region_start") or lead_pos
        r_end = row.get("region_end") or lead_pos

        # the window reaches window_size beyond whichever of lead and centroid is further out
        dist_lead_to_centroid = abs(centroid - lead_pos)
        search_radius = max(window_size, dist_lead_to_centroid + window_size)

        search_start = max(1, min(lead_pos, centroid) - search_radius)
        search_end = max(lead_pos, centroid) + search_radius

        # query Ensembl, skip the locus on failure
        try:
            genes_in_window = ensembl.genes_at_locus(
                contig=chrom, position=search_start, end=search_end
            )
        except Exception as e:
            logger.warning(
                f"Ensembl query failed for {locus_id} at {chrom}:{search_start}-{search_end}: {e}"
            )
            continue

        # keep each gene of the wanted biotypes once
        seen_genes = set()

        for gene in genes_in_window:
            if gene.biotype not in biotype or gene.gene_id in seen_genes:
                continue
            seen_genes.add(gene.gene_id)

            # TSS is the gene boundary on the transcribed strand
            tss_pos = gene.start if gene.strand == "+" else gene.end

            # --> lead distances <-- #
            # lead to TSS, unsigned and signed (negative = upstream)
            dist_lead_tss = abs(lead_pos - tss_pos)
            signed_dist_lead_tss = lead_pos - tss_pos

            # lead to footprint, 0 inside the gene body
            if gene.start <= lead_pos <= gene.end:
                dist_lead_footprint = 0
            else:
                dist_lead_footprint = min(
                    abs(gene.start - lead_pos), abs(gene.end - lead_pos)
                )

            # --> centroid distances <-- #
            # centroid to TSS, unsigned and signed
            dist_centroid_tss = abs(centroid - tss_pos)
            signed_dist_centroid_tss = centroid - tss_pos

            # centroid to footprint, 0 inside the gene body
            if gene.start <= centroid <= gene.end:
                dist_centroid_footprint = 0
            else:
                dist_centroid_footprint = min(
                    abs(gene.start - centroid), abs(gene.end - centroid)
                )

            # one row per (locus, gene)
            result_rows.append(
                {
                    "locus_id": locus_id,
                    "marker_name_hg38": marker_name,
                    "chr_hg38_ensembl": chrom,
                    "pos_hg38": lead_pos,
                    "ensembl_gene_id": gene.gene_id,
                    "hgnc_symbol": gene.gene_name,
                    "gene_biotype": gene.biotype,
                    "strand": gene.strand,
                    "TSS_position": tss_pos,
                    # region
                    "region_start": r_start,
                    "region_end": r_end,
                    "region_centroid_weighted": centroid,
                    # distances
                    "distance_to_footprint": dist_lead_footprint,
                    "distance_to_tss": dist_lead_tss,
                    "signed_distance_to_tss": signed_dist_lead_tss,
                    "weighted_distance_to_footprint": dist_centroid_footprint,
                    "weighted_distance_to_tss": dist_centroid_tss,
                    "signed_weighted_distance_to_tss": signed_dist_centroid_tss,
                }
            )

    # --> output <-- #
    if not result_rows:
        return None
    # sort by locus, then footprint distance, TSS distance breaks ties
    return (
        pl.DataFrame(result_rows)
        .sort(["locus_id", "distance_to_footprint", "distance_to_tss"])
    )
