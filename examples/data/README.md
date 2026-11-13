# Example input data

Input files for the two configs in `configs/example/`. They show the input formats
`run.py` accepts. Running them still needs the full reference data. The tests do not use
these files.

## Files

| File | Rows | What it shows |
| --- | --- | --- |
| `lead_variants_test.tsv` | 12 | lead variants, hg38, one row per locus |
| `proxies_test.tsv` | 62 | LD proxies linked to the leads by `locus_id`, named by `proxy_id` |
| `lead_variants_test_hg19.tsv` | 12 | the same leads in hg19 |
| `proxies_test_hg19.tsv` | 62 | the same proxies in hg19 |
| `phenotype_map_test.tsv` | 11 | locus -> phenotype index, for graph scoring |

