# LD reference panel

Stage 01, `pipelines/01_ld_proxies.py`, is a sample pipeline. It retrieves the LD proxies
of each lead from the 1000 Genomes EUR release on hg38 with plink2, and nothing else. It
is provided as is and is not maintained for other panels or builds. Another panel, such
as UK Biobank, or another build works only after manual changes to the script.

For a real analysis, build your own LD pipeline in line with your reference panel and
give its lead and proxy tables to stage 02 (section 1). The LD of the paper came from UK
Biobank through a separate pipeline and was given to DEFA this way.

## 1. Use your own lead and proxy tables

Set `lead_variant_file` and `proxy_file` for the study instead of `leads_file`. Stage 01
then skips the study and stage 02 reads your tables.

| File | Columns |
|---|---|
| lead | `locus_id`, `chr_hg38`, `pos_hg38`, `ref_hg38`, `alt_hg38` |
| proxy | `locus_id`, `proxy_id`, `r2`, `chr_hg38`, `pos_hg38`, `ref_hg38`, `alt_hg38` |

The `hg19` columns work too. Example rows and the format rules are in the Input section
of the [`README`](README.md#input), complete files in [`examples/data/`](examples/data/).

```yaml
studies:
  my_study:
    lead_variant_file: /path/to/leads_my_study.tsv
    proxy_file:        /path/to/proxies_my_study.tsv
```

## 2. Build the 1000 Genomes EUR panel (hg38)

Stage 01 reads the 525 unrelated EUR samples of 1000 Genomes. Download the hg38 files from
the [PLINK 2 resources page](https://www.cog-genomics.org/plink/2.0/resources):

- `all_hg38.pgen.zst`, `all_hg38.pvar.zst` and the matching `.psam`, saved as `all_hg38.psam`
- `deg2_hg38.king.cutoff.out.id`, the list of related samples to remove

```bash
cd /path/to/1kg_hg38

# plink2 reads .pvar.zst directly but needs the .pgen decompressed
plink2 --zst-decompress all_hg38.pgen.zst > all_hg38.pgen

# unrelated EUR samples, 525 of them
plink2 --pfile all_hg38 vzs \
       --remove deg2_hg38.king.cutoff.out.id \
       --keep-cat-pheno SuperPop --keep-cat-names EUR \
       --write-samples --out 1kg38_eur

# one fileset per chromosome: biallelic ACGT SNPs, EUR MAF >= 1%, no duplicate ids
mkdir -p 1kg38_eur_pfiles
for chr in {1..22} X; do
  plink2 --pfile all_hg38 vzs --chr $chr \
         --keep 1kg38_eur.id \
         --snps-only just-acgt --max-alleles 2 --maf 0.01 \
         --rm-dup exclude-all \
         --output-chr 26 \
         --make-pgen --out 1kg38_eur_pfiles/plink2_chr${chr}
done

# the decompressed .pgen is ~10 GB and only needed for the split
rm all_hg38.pgen
```

The finished panel is ~530 MB, built with plink2 v2.0.0-a.7.6LM.

```yaml
ld_proxies:
  panel:
    dir: /path/to/1kg_hg38/1kg38_eur_pfiles
```

The script assumes this exact panel:

- build hg38
- one fileset per chromosome, `plink2_chr{1..22,X}.{pgen,pvar,psam}`
- plain chromosome codes, `1` .. `22` and `X` or `23`
- all samples in the psam are used
- unphased r2 (`--r2-unphased`)

The plink2 binary, the window and the r2 cut-off are set in the `ld_proxies` block of
[`pipelines/config.template.yaml`](pipelines/config.template.yaml).

## 3. Give your leads

Stage 01 only needs the lead positions, in hg38:

```
locus_id     chrom_hg38  pos_hg38
rs515135     2           21063185
```

Stage 01 does not lift over, so lift your leads first. `locus_id` is a free label.

Leads are matched to the panel by chromosome and position, not by rsID or alleles.
Matching on the full variant id would be stricter and may replace this later. A lead
whose position is not in the panel is dropped, so check
`ld_proxies/leads_resolved_{study}.tsv` after each run:

| Status | Meaning | What to do |
|---|---|---|
| `matched` | one panel variant at the position | check it is your variant, an indel lead can match a SNP at the same position |
| `not_in_panel` | no panel variant at the position: an indel, EUR MAF below 1%, or no file for the chromosome (Y, MT) | use a tag variant the panel holds |
| `ambiguous_position` | more than one panel variant at the position | use a tag variant, the lead is skipped |
