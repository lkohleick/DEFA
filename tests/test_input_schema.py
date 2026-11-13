################################################
####   Unit tests: input files and loaders  ####
################################################
"""Tests for the id-column schema of the input files and the lead plus proxy loader."""

from defa import input_schema as schema
from defa import variant_data_processing as vdp


# both builds on every row, so process_variants_dataframe never lifts over
COORDS = "chr_hg38\tpos_hg38\tref_hg38\talt_hg38\tchr_hg19\tpos_hg19\tref_hg19\talt_hg19"
LEAD_1 = "chr1\t100\tA\tG\tchr1\t90\tA\tG"
LEAD_2 = "chr2\t200\tC\tT\tchr2\t190\tC\tT"
PROXY_1 = "chr1\t150\tG\tC\tchr1\t140\tG\tC"


def _write(path, lines):
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def _lead_and_proxy(tmp_path, lead_id="locus_id", proxy_ids="locus_id\tproxy_id"):
    lead = _write(tmp_path / "leads.tsv", [
        f"{lead_id}\t{COORDS}",
        f"1\t{LEAD_1}",
        f"L2\t{LEAD_2}",
    ])
    proxy = _write(tmp_path / "proxies.tsv", [
        f"{proxy_ids}\tr2\t{COORDS}",
        f"1\tP1\t0.8\t{PROXY_1}",
    ])
    return lead, proxy


################################################
####               the reader               ####
################################################


def test_old_proxy_names_are_renamed(tmp_path):
    _lead, proxy = _lead_and_proxy(tmp_path, proxy_ids="rsid_lead\trsid_proxy")

    df = schema.read_input_table(proxy, "proxy")

    assert {"locus_id", "proxy_id"} <= set(df.columns)
    assert not {"rsid_lead", "rsid_proxy"} & set(df.columns)


################################################
####         lead plus proxy input          ####
################################################


def test_leads_and_proxies_stack_into_one_row_per_variant(tmp_path):
    lead, proxy = _lead_and_proxy(tmp_path)

    df = vdp.load_lead_and_proxy(lead, proxy).sort(["locus_id", "variant_id"])

    assert df["locus_id"].to_list() == ["1", "1", "L2"]
    assert df["variant_id"].to_list() == ["1", "P1", "L2"]
    assert df["scaling_factor"].to_list() == [1.0, 0.8, 1.0]
    assert df["marker_name_hg38"].to_list() == ["chr1:100_A_G", "chr1:150_C_G", "chr2:200_C_T"]
