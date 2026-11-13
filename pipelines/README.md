# Workflow

```bash
cp pipelines/config.template.yaml my_run/config.yaml
# edit my_run/config.yaml

python pipelines/01_ld_proxies.py      --config my_run/config.yaml
python pipelines/02_defa_annotation.py --config my_run/config.yaml
python pipelines/03_gene_scoring.py    --config my_run/config.yaml
```

Activate the `defa` environment first, see [INSTALL.md](../INSTALL.md).

## The three stages

| Script | Input | Output |
|--------|-------|--------|
| `01_ld_proxies.py` | lead TSV per study (`locus_id`, `chrom_hg38`, `pos_hg38`) and the 1000 Genomes EUR hg38 panel | `leads_resolved_{study}.tsv`, `leads_{study}.tsv`, `proxies_{study}.tsv`, `variant_to_endpoint_map_{study}.tsv` |
| `02_defa_annotation.py` | those tables | `defa_variant_table_{study}.tsv`, every variant-gene pair with its features |
| `03_gene_scoring.py` | the annotated table | `scored_table_{study}.tsv`, `best_genes_{study}.tsv` |

Stage 01 is optional. A study that sets its own `lead_variant_file` plus `proxy_file`,
in the format 01 writes, is skipped by 01 and picked up by 02.

## Flags

Every script takes the same three:

| Flag | Meaning |
|------|---------|
| `--config PATH` | the workflow config, required |
| `--study NAME` | run one study instead of all of them |
| `--force` | redo a study whose outputs are already present |

Without `--force` a finished study is skipped, so a failed run is resumed by starting
the same command again.

## Configs

`configs/base_config.yaml` holds the defaults and the `/path/to` placeholders. Do not edit
it. `configs/base_config.local.yaml` is gitignored and is loaded on top of it. Start it
from `configs/base_config.local.template.yaml`.

```
configs/base_config.yaml                defaults, /path/to placeholders, do not edit
`-- configs/base_config.local.yaml      your machine paths, merged on top
    `-- one run layer, merged on top:
        |-- {output_dir}/{study}/defa/study_config.yaml   written by stage 02 from my_run/config.yaml
        `-- configs/example/test_config_hg38.yaml         run.py format, annotation only

templates, copied once:
  configs/base_config.local.template.yaml -> configs/base_config.local.yaml
  pipelines/config.template.yaml          -> my_run/config.yaml, read by stages 01-03
```

A layer may only set keys that `base_config.yaml` defines. `my_run/config.yaml` is not
merged itself: stage 02 takes its `annotation` block and the study entries and writes
them to `study_config.yaml`, stages 01 and 03 read their own blocks directly.

## Input format

The columns of each file are listed under Input in the [README](../README.md#input),
complete files are in [`examples/data/`](../examples/data/).

- All files are tab-separated with a header row.
- Chromosomes carry the `chr` prefix, `chr1`..`chr22`, `chrX`. Stage 01 input is the one
  place without it.
- `locus_id` in the proxy file must match the lead file.
- The lead does not have to be in the proxy file, it counts as its own variant at r2 1.
  A locus without proxies is fine.
- Proxies below `scaling_factor_threshold`, by default r2 0.5, are dropped.
- A row with an empty field is dropped.
- One build is enough, DEFA lifts over to the other. Use the `hg19` columns for hg19.
- Extra columns are ignored.
- `phenotype_idx` is the PrimeKG `node_index` from
  `defa_data/prime_kg_20220425/nodes.csv`, only used for graph scoring. Each study has
  one phenotype.

## Output layout

```
{output_dir}/
`-- {study}/
    |-- ld_proxies/
    |   |-- leads_resolved_{study}.tsv  # every input locus and why it matched or did not
    |   |-- leads_{study}.tsv
    |   |-- proxies_{study}.tsv
    |   |-- variant_to_endpoint_map_{study}.tsv
    |   `-- tmp/                        # SNP lists and raw plink2 output
    |-- defa/
    |   |-- study_config.yaml           # the DEFA config this study ran with
    |   |-- cache/                      # step caches and VEP temp files
    |   |-- defa_{study}.log
    |   `-- defa_variant_table_{study}.tsv
    `-- scoring/
        |-- scored_table_{study}.tsv         # every variant-gene pair with its model scores and DEFA-mean
        `-- best_genes_{study}.tsv          # the top-ranked gene per locus, per model and DEFA-mean
```

Each model score is min-max scaled within its locus, `rel_<model> = (p - min) / (max - min)`,
so the top gene scores 1 and the lowest 0, and a locus with one gene scores 1. DEFA-mean,
`defa_mean_with_graph_features` and `defa_mean_without_graph_features`, is the mean of the
five domain-informed `rel_` scores of one family, and reaches 1 only when all five models
agree.

## Runtime

Stage 01 runs plink2 once per chromosome and takes minutes for a few hundred loci.
Stage 02 is the long one: it queries VEP, eQTL catalogues and the regulatory datasets,
and its memory grows with the number of loci. Stage 03 is inference only and
takes seconds.

On a cluster, submit each script as a single job, for example

```bash
sbatch --cpus-per-task=16 --mem=128G --time=24:00:00 \
       --wrap "python pipelines/02_defa_annotation.py --config my_run/config.yaml"
```

## Annotation only

`run.py` runs only the annotation stage on one dataset and writes
`defa_variant_table_{name}.tsv`. It does not find proxies or score genes. It takes a config
in the format of `configs/example/`:

```bash
python run.py --config configs/example/test_config_hg38.yaml
python run.py --config my_config.yaml --cache-dir cache/run_42   # one cache per parallel job
```

## Testing

`pipelines/test/` has a small test study of three coronary artery disease loci. How to
run it is in section 7 of [INSTALL.md](../INSTALL.md).
