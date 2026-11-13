<p align="center">
  <img src="DEFA.png" alt="DEFA" width="140">
</p>

<h1 align="center">DEFA</h1>

<p align="center">
  Domain-informed EFfector gene Assignment
</p>

<p align="center">
  <a href="https://github.com/lkohleick/DEFA/actions/workflows/test.yml">
    <img src="https://github.com/lkohleick/DEFA/actions/workflows/test.yml/badge.svg" alt="tests">
  </a>
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-blue.svg" alt="MIT license">
  </a>
  <img src="https://img.shields.io/badge/python-3.10%20%7C%203.11-blue.svg" alt="Python 3.10 | 3.11">
</p>

<p align="center">
  <a href="INSTALL.md"><b>Installation</b></a> |
  <a href="REFERENCE_DATA.md"><b>Reference data</b></a> |
  <a href="pipelines/README.md"><b>Workflow</b></a>
</p>

DEFA (Domain-informed EFfector gene Assignment) assigns effector genes to GWAS lead
variants. It annotates every candidate gene around each lead variant and its LD proxies
with functional genomic evidence and, when the trait is given, with its connectivity to
the trait in a biomedical knowledge graph. Several XGBoost classifiers, each trained on
quasi-true positives from a different line of evidence, then score the candidates. Their
scores are combined, and each one is also reported, so the assumption behind every
assignment is made explicit.

```
lead variants --01--> LD proxies --02--> annotated variant-gene table --03--> ranked genes
```

## What you need

| | |
|---|---|
| OS | Linux. Not tested on macOS or Windows |
| Software | micromamba, mamba or conda, two environments (`defa`, `vep`), plink2 only for stage 01 |
| Disk | ~125 GB: VEP cache + CADD ~105 GB, reference data ~19 GB, pyensembl ~1 GB |
| CPU | as many cores as set in `vep.fork`, default 8 |
| Memory | a few GB for a few hundred loci, over 100 GB for ~1000 loci |

## Quickstart

Install the software and reference data and run the three-locus test study, sections 1-7
of [`INSTALL.md`](INSTALL.md). The reference data and trained models are on Zenodo,
[doi:10.5281/zenodo.24160786](https://doi.org/10.5281/zenodo.24160786).
Then, from the repository root with `defa` active:

```bash
mkdir -p my_run
cp pipelines/config.template.yaml my_run/config.yaml
# set output_dir, one entry under studies:, and ld_proxies.panel.dir

python pipelines/01_ld_proxies.py      --config my_run/config.yaml
python pipelines/02_defa_annotation.py --config my_run/config.yaml
python pipelines/03_gene_scoring.py    --config my_run/config.yaml
```

With leads and proxies from your own LD pipeline, skip stage 01, see [Input](#input).
Flags, resuming, output layout and cluster resources are in
[`pipelines/README.md`](pipelines/README.md).

## Configs

| File | Holds | Who edits it |
|---|---|---|
| `configs/base_config.local.yaml` | where the reference data and VEP live on this machine | you, once per machine |
| `my_run/config.yaml` | studies, input files, LD panel, annotation and scoring options | you, once per analysis |

Do not edit `configs/base_config.yaml`, it holds the defaults. How the layers merge is in
[`pipelines/README.md`](pipelines/README.md#configs).

## Input

Stage 01 is only an example LD pipeline to get the workflow running. It retrieves LD
proxies from the 1000 Genomes EUR hg38 release, see [`LD_PANEL.md`](LD_PANEL.md). For a
real analysis, amend it or replace it with a more complex LD backbone that fits your
study. It only needs the lead positions, hg38, chromosomes `1`..`22`, `X`:

```
locus_id     chrom_hg38  pos_hg38
rs515135     2           21063185
```

If you already have leads and proxies, skip stage 01. Give the study a `lead_variant_file`
and `proxy_file` instead of `leads_file`, drop the `ld_proxies` block and run only stages
02 and 03. For graph scoring also set `gwas_endpoint_file`. All files are tab-separated
with a header row, chromosomes carry the `chr` prefix:

| File | Columns |
|---|---|
| lead | `locus_id`, `chr_hg38`, `pos_hg38`, `ref_hg38`, `alt_hg38` |
| proxy | `locus_id`, `proxy_id`, `r2`, `chr_hg38`, `pos_hg38`, `ref_hg38`, `alt_hg38` |
| endpoint map (graph scoring only) | `locus_id`, `phenotype_idx` |
| true-positive links (optional) | `locus_id`, `ensembl_gene_id` |

`phenotype_idx` is the PrimeKG `node_index` of the study phenotype. Look it up in
`defa_data/prime_kg_20220425/nodes.csv`, which ships with the DEFA reference data, e.g.
36187 for coronary artery disease. For several phenotypes on the same leads, set up one
study per phenotype. Complete files are in [`examples/data/`](examples/data/), the full
format rules in [`pipelines/README.md`](pipelines/README.md#input-format).

## Output

Per study, under `{output_dir}/{study}/`:

| File | Holds |
|---|---|
| `scoring/best_genes_{study}.tsv` | the top-ranked gene per locus, for each model and DEFA-mean |
| `scoring/scored_table_{study}.tsv` | every variant-gene pair with its model scores and DEFA-mean |
| `defa/defa_variant_table_{study}.tsv` | every variant-gene pair with all annotation columns |
| `ld_proxies/` | resolved leads, proxies, endpoint map |

Each model score is min-max scaled within its locus, so the top gene scores 1. DEFA-mean
is the mean of the five domain-informed scores of one family.
[`REFERENCE_DATA.md`](REFERENCE_DATA.md) describes each annotation column and its source.

## License

MIT, see [LICENSE](LICENSE). The license covers the DEFA code only. The reference data and
the external tools DEFA runs keep their own licenses, listed in
[`REFERENCE_DATA.md`](REFERENCE_DATA.md).
