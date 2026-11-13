# Reference data and annotation columns

All datasets with versions and sources are listed in Supplementary Table 2 of the
[preprint](<preprint DOI>).

### Gene and distance

| Column(s) | Dataset | Config key | Where |
|---|---|---|---|
| `distance_to_footprint`, `distance_to_tss`, `weighted_distance_to_footprint`, `weighted_distance_to_tss` | computed on Ensembl 111 gene models | `parameters.ensembl_release` | pyensembl cache |
| `hgnc_symbol`, `ensembl_gene_id` | Ensembl 111 | `parameters.ensembl_release` | pyensembl cache |

### Variant effect

| Column(s) | Dataset | Config key | Where |
|---|---|---|---|
| `VEP_max`, `VEP_sum`, `VEP_weighted_mean` | Ensembl VEP 106, SpliceRegion plugin | `vep.executable` | VEP install |
| `CADD_phred_max/mean/sum`, `CADD_raw_max/mean/sum`, `CADD_phred_weighted_mean`, `CADD_raw_weighted_mean` | CADD v1.6, gnomAD r3.0 indels | `vep.cache_dir` | user download |

### eQTL

| Column(s) | Dataset | Config key | Where |
|---|---|---|---|
| `eqtl_<tissue>_scaled_sum`, `eqtl_<tissue>_weighted_mean` (99 tissues and cell types) | eQTL Catalogue release 7, max SuSiE PIP, p <= 5e-8, scaled by the variant `scaling_factor` and summed or weighted-mean per locus | `paths.eqtl.eqtl_catalogue_max_pip` | DEFA Zenodo |
| `eqtl_weight_sum` | eQTL Catalogue release 7, summed `scaling_factor` of the locus variants with eQTL evidence for the gene | `paths.eqtl.eqtl_catalogue_max_pip` | DEFA Zenodo |
| `CLPP_GTEX_eQTL_sum/max` | GTEx v8 CAVIAR fine-mapping | `paths.eqtl.gtex_caviar_flames` | FLAMES bundle |
| `CLPP_eQTLgen_sum/max` | eQTLGen whole blood | `paths.eqtl.eqtlgen` | FLAMES bundle |
| `max_CLPP_eQTL_eQTLCatalog_*`, `mean_CLPP_eQTL_eQTLCatalog_*`, rQTL equivalents | eQTL Catalogue eQTL and rQTL | `paths.eqtl.eqtl_catalogue_eqtl`, `paths.eqtl.eqtl_catalogue_rqtl` | FLAMES bundle |

### Regulatory and chromatin

| Column(s) | Dataset | Config key | Where |
|---|---|---|---|
| `ABC_score` | CATlas v1 cCREs and ABC scores, Ensembl 112 transcript -> gene map | `paths.gene_regulation.catlas.*`, `paths.gene_information.enst_to_ensg_map` | DEFA Zenodo |
| `ABC_EP_sum/max/mean` | ABC enhancer-promoter links | `paths.gene_regulation.abc_ep_data` | FLAMES bundle |
| `ABC_CRISPR_sum/max/mean/weighted_mean` | ABC CRISPRi-validated links | `paths.gene_regulation.abc_crispr_data` | FLAMES bundle |
| `HACER_CAGE_sum/max/mean`, `HACER_PRO-seq_GRO-seq_sum/max/mean` | HACER enhancer catalogue | `paths.gene_regulation.hacer` | FLAMES bundle |
| `Jung_HiC_sum/max/mean` | Jung et al. promoter capture Hi-C | `paths.gene_regulation.hic_annotation.jung_hic_data` | FLAMES bundle |
| `Javierre_HiC_sum/max/mean` | Javierre et al. promoter capture Hi-C | `paths.gene_regulation.hic_annotation.javierre_hic_data` | FLAMES bundle |
| `F5_ep_sum/max` | FANTOM5 enhancer-promoter correlation | `paths.gene_regulation.fantom5` | FLAMES bundle |
| `GeneHancer_sum/max` | GeneHancer | `paths.gene_regulation.genehancer` | FLAMES bundle |
| `RoadmapEpi_sum/max` | Roadmap Epigenomics | `paths.gene_regulation.roadmap` | FLAMES bundle |
| `EpiMap_sum/max` | EpiMap | `paths.gene_regulation.epimap` | FLAMES bundle |
| `cicero_sum/max/mean` | Cicero whole-blood co-accessibility | `paths.gene_regulation.cicero` | FLAMES bundle |
| `Promoter_sum/max` | promoter regions per transcript | `paths.gene_regulation.promoter` | FLAMES bundle |

### Disease and drug targets

| Column(s) | Dataset | Config key | Where |
|---|---|---|---|
| `omim_gene` (0/1) | OMIM `mim2gene.txt` | `paths.omim.genes_catalogue` | user download |
| `open_targets_advanced_drug_target` (0/1) | Open Targets Platform, 2024-07-02 | `paths.open_targets.drug_targets` | DEFA Zenodo |

### Graph scores

| Column(s) | Dataset | Config key | Where |
|---|---|---|---|
| `graph_distance`, `adamic_adar`, `pagerank_score` (+ `*_with_exclusions`) | PrimeKG 2022-04-25 | `graph_scorer.primekg_dir`, `graph_scorer.ensembl_to_primekg_idx_file` | DEFA Zenodo |
| `txgnn_distmult` | TxGNN DistMult embeddings | `graph_scorer.txgnn_export_dir` | DEFA Zenodo |
