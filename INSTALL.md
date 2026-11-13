# Installation

To run DEFA you need two conda environments, VEP with its cache and the reference
datasets. plink2 and an LD panel are only needed if you retrieve LD proxies with stage 01.
This guide goes through them in order and ends with a small test study. What each
dataset contains, and its license, is in [`REFERENCE_DATA.md`](REFERENCE_DATA.md).

The commands use micromamba, but mamba or conda work the same way.

## What you need

| # | Component | Size | Needed for |
|---|---|---|---|
| 1 | `defa` environment | | everything |
| 2 | pyensembl cache | ~1 GB | gene annotation |
| 3 | `vep` environment | | every run |
| 4 | VEP cache, FASTA and CADD scores | ~105 GB | every run |
| 5 | `plink2` | | optional, only to retrieve proxies with `pipelines/01_ld_proxies.py` |
| 6 | reference datasets | 3.7 GB download, ~19 GB unpacked | annotation and scoring |
| 7 | 1000 Genomes EUR panel | ~530 MB | optional, only to retrieve proxies with `pipelines/01_ld_proxies.py` |


## 1. The `defa` environment

```bash
git clone https://github.com/lkohleick/DEFA.git
cd DEFA

micromamba env create -f environment.yml
micromamba activate defa
```

Then download the Ensembl gene annotation. The release has to be 111, the same as
`parameters.ensembl_release` in `configs/base_config.yaml`:

```bash
pyensembl install --release 111 --species homo_sapiens
```

This goes to `~/.cache/pyensembl`. If your home directory has a quota, set
`PYENSEMBL_CACHE_DIR` to somewhere else first.

## 2. VEP

### 2a. The `vep` environment

```bash
micromamba create -n vep -c conda-forge -c bioconda ensembl-vep=106
micromamba env list          # note the path, it goes in the config
```

The experiments used VEP 106. You can try a newer release, then use its number in every
place that says 106: `ensembl-vep=`, `--CACHE_VERSION` and the `106_GRCh38` folder in
`vep.fasta_file`.

### 2b. Cache and reference FASTA

Choose one directory for all VEP data. This becomes `vep.cache_dir` in the config.

```bash
export VEP_DIR=/path/to/vep_cache
mkdir -p "$VEP_DIR"
micromamba activate vep

vep_install -a cf -s homo_sapiens -y GRCh38 \
            --CACHE_VERSION 106 \
            -c "$VEP_DIR" \
            --CONVERT
```

This is about 25 GB and creates `$VEP_DIR/homo_sapiens/106_GRCh38/`.

### 2c. Plugins

DEFA uses the CADD and SpliceRegion plugins. They have to go into `.vep/Plugins38`
inside the cache directory:

```bash
vep_install -a p --PLUGINS CADD,SpliceRegion \
            --PLUGINSDIR "$VEP_DIR/.vep/Plugins38"
```

### 2d. CADD scores

The CADD plugin needs the v1.6 GRCh38 score files, about 80 GB in total. Download them
into this exact path. If the link no longer works, check
https://cadd.gs.washington.edu/download.

```bash
CADD_DIR="$VEP_DIR/vep_plugins/CADD/krishna.gs.washington.edu/download/CADD/v1.6/GRCh38"
mkdir -p "$CADD_DIR"
cd "$CADD_DIR"

BASE=https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38
for f in whole_genome_SNVs.tsv.gz whole_genome_SNVs.tsv.gz.tbi \
         gnomad.genomes.r3.0.indel.tsv.gz gnomad.genomes.r3.0.indel.tsv.gz.tbi; do
    wget -c "$BASE/$f"
done
```

Make sure the `.tbi` index files are there too, otherwise VEP fails when it runs. `wget -c`
picks up where it stopped if the download is interrupted.

CADD offers the same files from two servers, the University of Washington in Seattle (used
above) and the Berlin Institute of Health at Charite. From Europe the Berlin server is much
faster:

```bash
BASE=https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38
```

Keep `CADD_DIR` as it is either way, VEP expects that path.

### 2e. Final layout

Everything has to be in this one directory:

```
{vep.cache_dir}/
|-- homo_sapiens/106_GRCh38/
|   `-- Homo_sapiens.GRCh38.dna.toplevel.fa.gz
|-- .vep/Plugins38/
`-- vep_plugins/CADD/krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/
    |-- whole_genome_SNVs.tsv.gz
    |-- whole_genome_SNVs.tsv.gz.tbi
    |-- gnomad.genomes.r3.0.indel.tsv.gz
    `-- gnomad.genomes.r3.0.indel.tsv.gz.tbi
```

VEP only runs on GRCh38, so leave `parameters.genome_build` at `hg38`. Input files in hg19
are fine, DEFA lifts them over to hg38 when it loads them. Set `vep.fork` to the number of
cores you have (default 8).

## 3. plink2 (optional)

Only `pipelines/01_ld_proxies.py` needs plink2, to retrieve LD proxies for your leads.
Skip this section if you give your own lead and proxy tables. The test study in section
7 ships its proxies and does not need it. Download it from
https://www.cog-genomics.org/plink/2.0/ or install it into the `defa` environment:

```bash
micromamba install -c bioconda plink2
```

If it is not on your `PATH`, set `ld_proxies.plink2_binary` in the workflow config.

## 4. Reference datasets

There are three downloads:

| Directory | Source | Download | Config key |
|---|---|---|---|
| `defa_data/` | DEFA, [doi:10.5281/zenodo.24160786](https://doi.org/10.5281/zenodo.24160786) | `defa_data_20261112.zip`, 2.0 GB | `defa_data_dir` |
| `Annotation_data/` | FLAMES, [doi:10.5281/zenodo.12635505](https://doi.org/10.5281/zenodo.12635505) | `Annotation_data.tar.gz` only, 1.7 GB | `flames_data_dir` |
| `omim/` | OMIM, https://omim.org/downloads | `mim2gene.txt` | `paths.omim.genes_catalogue` |

At the end it should look like this:

```
$REF/
|-- defa_data/
|   |-- CATlas/
|   |-- ensembl/
|   |-- eqtl_catalogue/
|   |-- open_targets/
|   |-- models/defa_models_v1.0.0/
|   |-- prime_kg_20220425/
|   |-- primekg/
|   `-- txgnn_export_numpy/
|-- Annotation_data/
`-- omim/mim2gene.txt
```

### 4a. DEFA data

Download the data from https://doi.org/10.5281/zenodo.24160786 and unpack it:

```bash
REF=/path/to/reference_data
mkdir -p "$REF"

unzip /path/to/defa_data_20261112.zip -d /path/to/defa_download
cd /path/to/defa_download
for f in defa_*.tar*; do tar -xf "$f" -C "$REF"; done
```

### 4b. FLAMES

You only need `Annotation_data.tar.gz`:

```bash
cd /path/to/flames_download
tar -xzf Annotation_data.tar.gz -C "$REF"
```


### 4c. OMIM

```bash
mkdir -p "$REF/omim"
cp /path/to/mim2gene.txt "$REF/omim/"
```

### 4d. Config

Your paths go in `configs/base_config.local.yaml`. It is gitignored and loaded on top of
`configs/base_config.yaml`. Copy the template that comes with the repository:

```bash
cp configs/base_config.local.template.yaml configs/base_config.local.yaml
```

Open `configs/base_config.local.yaml` and replace every `/path/to/...` with your own
path: the reference data folder from 4a, the VEP cache from 2b and the `vep` environment
from 2a. Any path still starting with `/path/to` is reported by the check in section 6.

All other reference paths are found through `defa_data_dir` and `flames_data_dir`. A key
that `base_config.yaml` does not define is rejected.

## 5. LD reference panel (optional)

`pipelines/01_ld_proxies.py` is a sample pipeline that retrieves LD proxies from the 1000
Genomes EUR release on hg38. It is provided as is and works with other panels or builds
only after manual changes to the script. How to build the panel is in section 2 of
[`LD_PANEL.md`](LD_PANEL.md). It is open access, ~530 MB. The test study ships proxies
taken from it.

The panel has 525 EUR samples and is there to get the workflow running. For a real
analysis, build your own LD pipeline in line with your reference panel, skip stage 01
and give its lead and proxy tables, section 1 of [`LD_PANEL.md`](LD_PANEL.md).

```yaml
ld_proxies:
  panel:
    dir: /path/to/1kg_hg38/1kg38_eur_pfiles
```

Leads are matched to the panel by chromosome and position, not by rsID. Section 3 of
[`LD_PANEL.md`](LD_PANEL.md) explains how to give your leads so none get lost.

## 6. Check the install

```bash
python run.py --config configs/example/test_config_hg38.yaml --check-only
```

This checks that all reference paths, the VEP paths, the plugins and the CADD files can
be found. Fix whatever it reports before running anything.

To run the unit tests:

```bash
pytest tests/
```

## 7. First run: the test study

`pipelines/test/` has a small coronary artery disease study with three non-coding leads
where the causal gene is well known: rs515135 (APOB), rs12149545 (CETP) and rs894211
(LPL). The phenotype is the PrimeKG node for coronary artery disease, 36187, so graph
scoring can be turned on too. The leads, their EUR LD proxies with r2 >= 0.5 and the
endpoint map ship in `pipelines/test/`, so stage 01 and the LD panel are not needed. Stages
02 and 03 run on the real reference data, so sections 1, 2 and 4 have to be set up first.

Graph scoring is on in the test config. On the first run DEFA builds the PrimeKG graph
with NetworkX, which takes a while, and caches it at `graph_scorer.cache_file`, by
default `cache/graph_scorer_graph.cache`. Later runs load the cache. To run without
graph features, set `run_graph_scorer: false` under `annotation:`.

Then, from the repository root with `defa` active:

```bash
python pipelines/02_defa_annotation.py --config pipelines/test/config.yaml
python pipelines/03_gene_scoring.py    --config pipelines/test/config.yaml
```

The output goes to `pipelines/test/output/test_run/`:

| File | Written by | Holds |
|---|---|---|
| `defa/defa_variant_table_test_run.tsv` | 02 | all variant-gene pairs of the three loci with their annotations |
| `scoring/best_genes_test_run.tsv` | 03 | the top gene per locus |

In `best_genes_test_run.tsv` almost all models should rank APOB, CETP and LPL first.
The whole run takes about 2.5 minutes with graph scoring on, plus the graph build on
the first run.

If a step fails, rerun the same command and it continues from where it stopped. Use
`--force` to redo a finished step. To delete the output, run `rm -r pipelines/test/output`.
