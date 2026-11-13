################################################
####        OMIM: disease gene flag         ####
################################################
"""Flag candidate genes listed in the OMIM gene catalogue."""

import polars as pl


def create_omim_gene_annotation(
    genes_in_region: pl.DataFrame, omim_genes_catalogue_path: str
) -> pl.DataFrame:
    """`omim_gene` (0/1) per candidate gene, from OMIM phenotype and gene entries with an Ensembl id."""
    # --> read the catalogue <-- #
    dtypes = {
        "MIM Number": pl.Int32,
        "MIM Entry Type (see FAQ 1.3 at https://omim.org/help/faq)": pl.Categorical,
        "Entrez Gene ID (NCBI)": pl.Int32,
        "Approved Gene Symbol (HGNC)": pl.Utf8,
        "Ensembl Gene ID (Ensembl)": pl.Utf8,
    }

    # the first 5 rows are comments, unparseable values become null
    omim_genes_catalogue = pl.read_csv(
        omim_genes_catalogue_path,
        separator="\t",
        skip_rows=5,
        has_header=True,
        dtypes=dtypes,
        null_values=["", "NA", "None"],
        ignore_errors=True,
    )

    # --> keep phenotype and gene entries with an Ensembl id <-- #
    valid_entry_types = [
        "predominantly phenotypes",
        "phenotype",
        "gene/phenotype",
        "gene",
    ]
    omim_genes_catalogue = omim_genes_catalogue.filter(
        (
            pl.col("MIM Entry Type (see FAQ 1.3 at https://omim.org/help/faq)").is_in(
                valid_entry_types
            )
        )
        & (pl.col("Ensembl Gene ID (Ensembl)").str.len_chars() > 0)
    )

    # --> flag candidate genes <-- #
    ensembl_ids_in_omim = omim_genes_catalogue.select(
        "Ensembl Gene ID (Ensembl)"
    ).to_series()

    omim_gene_annotation = genes_in_region.with_columns(
        [pl.col("ensembl_gene_id").is_in(ensembl_ids_in_omim).cast(pl.Int8).alias("omim_gene")]
    )
    omim_gene_annotation = omim_gene_annotation.select(
        [
            "ensembl_gene_id",
            "hgnc_symbol",
            "omim_gene",
        ]
    )
    return omim_gene_annotation
