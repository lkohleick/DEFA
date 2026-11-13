################################################
####    Variant table: join all evidence    ####
################################################
"""Join every annotation source onto the candidate genes and emit the variant table."""

import polars as pl

from .input_schema import read_input_table


################################################
####            evidence tables             ####
################################################


def define_eqtl_evidence_for_table(eqtl_catalogue_annotation, variant_data):
    """eQTL catalogue tissue columns, rekeyed from locus_id to marker_name_hg38."""
    eqtl_catalogue_annotation = eqtl_catalogue_annotation.join(
        variant_data["locus_id", "marker_name_hg38"],
        left_on="locus_id",
        right_on="locus_id",
        how="left",
    )
    eqtl_catalogue_annotation = eqtl_catalogue_annotation.drop("locus_id")
    # prefix tissue columns with eqtl_ and squash spaces, they end up as column names
    eqtl_catalogue_annotation = eqtl_catalogue_annotation.rename(
        {
            col: f"eqtl_{col.replace(' ', '_')}"
            for col in eqtl_catalogue_annotation.columns
            if col not in ["marker_name_hg38", "gene_id"]
        }
    )

    return eqtl_catalogue_annotation


def define_catlas_evidence_for_table(
    catlas_annotation,
    variant_data,
    threshold_abc=0.1,
):
    """cATLAS ABC scores above threshold_abc, one row per (marker, gene)."""
    catlas_annotation = catlas_annotation.select(["locus_id", "gene_id", "ABC_score"])

    # Filter by ABC score threshold
    catlas_annotation = catlas_annotation.filter(pl.col("ABC_score") >= threshold_abc)

    # one score per (locus, gene): keep="first" takes the first row, not the strongest
    catlas_annotation = catlas_annotation.unique(
        subset=["locus_id", "gene_id"], keep="first"
    )

    # rekey to marker name
    catlas_annotation = catlas_annotation.join(
        variant_data[["locus_id", "marker_name_hg38"]], on="locus_id", how="left"
    )

    # drop the loci that carry no marker name
    catlas_annotation = catlas_annotation.drop("locus_id")
    catlas_annotation = catlas_annotation.filter(
        pl.col("marker_name_hg38").is_not_null()
    )

    return catlas_annotation.select(["marker_name_hg38", "gene_id", "ABC_score"])


def define_true_positive_column(config):
    """True-positive variant-gene links, flagged with a constant true_positive of 1."""
    true_positive_links = read_input_table(
        config.dataset.true_positive_links_file, "true_positive"
    )
    # the file lists only true links, so the flag is constant; 0s come from fill_null
    return true_positive_links.with_columns(pl.lit(1).alias("true_positive"))


################################################
####          assemble the table            ####
################################################


def create_variant_table(
    variant_data,
    vep_annotations,
    genes_in_region,
    eqtl_catalogue_annotation,
    catlas_annotation,
    ensembl_gene_info_df,
    omim_genes_annotation,
    ot_drug_targets_annotation,
    flames_annotations,
    config,
):
    """One row per (marker, gene) carrying every annotation block, sorted by locus."""
    # skeleton: candidate genes plus the four distance features
    variant_table = genes_in_region.select(
        [
            "locus_id",
            "marker_name_hg38",
            "ensembl_gene_id",
            "hgnc_symbol",
            "distance_to_footprint",
            "distance_to_tss",
            "weighted_distance_to_footprint",
            "weighted_distance_to_tss",
        ]
    )
    # the only dedupe on locus_id, everything below keys on marker_name_hg38
    variant_table = variant_table.unique(subset=["locus_id", "ensembl_gene_id"])

    # add VEP Evidence
    variant_table = variant_table.join(
        vep_annotations[
            [
                "locus_id",
                "ensembl_gene_id",
                "VEP_max",
                "VEP_sum",
                "CADD_phred_max",
                "CADD_phred_mean",
                "CADD_phred_sum",
                "CADD_raw_max",
                "CADD_raw_mean",
                "CADD_raw_sum",
                "VEP_weighted_mean",
                "CADD_phred_weighted_mean",
                "CADD_raw_weighted_mean",
            ]
        ],
        left_on=["locus_id", "ensembl_gene_id"],
        right_on=["locus_id", "ensembl_gene_id"],
        how="left",
    )
    # --> rekey both locus-level tables to marker_name_hg38 <-- #
    eqtl_catalogue_annotation = define_eqtl_evidence_for_table(
        eqtl_catalogue_annotation,
        variant_data,
    )
    catlas_annotation = define_catlas_evidence_for_table(
        catlas_annotation,
        variant_data,
    )
    # OMIM and Open Targets are gene-level, they join on the gene alone
    omim_genes_annotation = omim_genes_annotation.select(
        ["ensembl_gene_id", "omim_gene"]
    )
    ot_drug_targets_annotation = ot_drug_targets_annotation.select(
        ["ensembl_gene_id", "open_targets_advanced_drug_target"]
    )
    # add true-positive variant-gene links
    if config.dataset.true_positive_links_file:
        true_positive_links = define_true_positive_column(config)
        variant_table = variant_table.join(
            true_positive_links["locus_id", "ensembl_gene_id", "true_positive"],
            on=["locus_id", "ensembl_gene_id"],
            how="left",
        )
        variant_table = variant_table.unique(
            subset=["marker_name_hg38", "ensembl_gene_id"]
        )

        # the two counts must agree, a shortfall means the dedupe dropped a row
        print(
            "Number of original true links", true_positive_links["true_positive"].sum()
        )
        print("Number of true links in variant table", variant_table["true_positive"].sum())
    else:
        pass

    # add catlas annotation
    variant_table = variant_table.join(
        catlas_annotation,
        left_on=["marker_name_hg38", "ensembl_gene_id"],
        right_on=["marker_name_hg38", "gene_id"],
        how="left",
    )
    variant_table = variant_table.unique(subset=["marker_name_hg38", "ensembl_gene_id"])

    # add flames annotation
    for annotation_name, annotation_data in flames_annotations.items():
        if annotation_data.is_empty():
            print(f"Adding zero-filled columns for empty annotation: {annotation_name}")

            # an empty frame still carries its schema, so take the names off it
            key_columns = {"locus_id", "ensembl_gene_id"}
            annotation_columns = [
                col for col in annotation_data.columns if col not in key_columns
            ]

            # write the zeros directly, a join against no rows would drop the columns
            for col in annotation_columns:
                variant_table = variant_table.with_columns(
                    pl.lit(0.0).alias(col).cast(pl.Float64)
                )
        else:
            print(f"Adding annotation data for: {annotation_name}")
            variant_table = variant_table.join(
                annotation_data,
                left_on=["locus_id", "ensembl_gene_id"],
                right_on=["locus_id", "ensembl_gene_id"],
                how="left",
            )

            # genes the annotation never saw score 0, not null
            key_columns = {"locus_id", "ensembl_gene_id"}
            annotation_columns = [
                col for col in annotation_data.columns if col not in key_columns
            ]

            for col in annotation_columns:
                variant_table = variant_table.with_columns(pl.col(col).fill_null(0.0))

    variant_table = variant_table.unique(subset=["marker_name_hg38", "ensembl_gene_id"])

    # add OMIM annotations
    variant_table = variant_table.join(
        omim_genes_annotation, on="ensembl_gene_id", how="left"
    )
    variant_table = variant_table.unique(subset=["marker_name_hg38", "ensembl_gene_id"])

    # add OT Drug Target annotation
    variant_table = variant_table.join(
        ot_drug_targets_annotation, on="ensembl_gene_id", how="left"
    )
    variant_table = variant_table.unique(subset=["marker_name_hg38", "ensembl_gene_id"])

    # add eQTL catalogue tissue evidence
    variant_table = variant_table.join(
        eqtl_catalogue_annotation,
        left_on=["marker_name_hg38", "ensembl_gene_id"],
        right_on=["marker_name_hg38", "gene_id"],
        how="left",
    )
    variant_table = variant_table.unique(subset=["marker_name_hg38", "ensembl_gene_id"])

    # blanket fill: this is where true_positive picks up its 0s
    variant_table = variant_table.fill_null(0)
    variant_table = (
        variant_table.unique()
    )

    # sort variant table after locus id and distance to lead
    variant_table = variant_table.sort(
        ["locus_id", "distance_to_footprint"], descending=[False, False]
    )

    return variant_table
