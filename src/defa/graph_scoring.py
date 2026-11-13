################################################
####    Graph scoring: PrimeKG and TxGNN    ####
################################################
"""PrimeKG phenotype-to-gene graph scores and TxGNN DistMult scores, joined onto the variant table."""

import os
import pickle
import networkx as nx
import pandas as pd
import polars as pl
from tqdm import tqdm
from loguru import logger
import numpy as np
import json


################################################
####            txgnn embeddings            ####
################################################


class TxGNNOfflineScorer:
    """DistMult scores from an offline TxGNN export of node embeddings and relation weights."""
    def __init__(self, export_dir):
        self.export_dir = export_dir

        # relation weights [R, D] and their index
        self.W = np.load(os.path.join(export_dir, "relations", "W.npy"))
        self.rel2idx = json.load(open(os.path.join(export_dir, "relations", "rel2idx.json")))

        # per node type: embeddings [N, D] and id -> row index
        self._emb = {}
        self._id2idx = {}

    def _key(self, ntype):
        """Node type as a directory name."""
        return ntype.replace("/", "_")

    def load_type(self, ntype):
        """Load embeddings and the id index of one node type, once."""
        if ntype in self._emb:
            return
        base = os.path.join(self.export_dir, "nodes", self._key(ntype))
        emb = np.load(os.path.join(base, "emb.npy"))
        idx2id = json.load(open(os.path.join(base, "idx2id.json")))
        id2idx = {v: int(k) for k, v in idx2id.items()}
        self._emb[ntype] = emb
        self._id2idx[ntype] = id2idx

    def distmult(self, src_type, src_id, rel, dst_type, dst_id):
        """sum(h * r * t) for one pair, None when either node has no embedding."""
        self.load_type(src_type)
        self.load_type(dst_type)

        # row indices of both nodes
        si = self._id2idx[src_type].get(str(src_id))
        ti = self._id2idx[dst_type].get(str(dst_id))
        if si is None or ti is None:
            return None

        # relation, head and tail vectors, each [D]
        r = self.W[self.rel2idx[rel]]
        h = self._emb[src_type][si]
        t = self._emb[dst_type][ti]
        return float(np.sum(h * r * t))


################################################
####              graph scorer              ####
################################################


class GraphScorer:
    """PrimeKG as an undirected NetworkX graph, scored from phenotype nodes to gene nodes."""
    def __init__(self):
        # graph, gene node set, and node sets per phenotype type
        self.G = None
        self.gene_set = None
        self.phenotype_sets = None

    # --> graph and cache <-- #
    def build_graph(self, primekg_df):
        """Build the NetworkX graph and the node sets from the PrimeKG edge table."""
        print("Building graph...")
        G = nx.Graph()

        def normalize_id(v):
            """Node id as a string, with a float-style `.0` suffix dropped."""
            if v is None or (isinstance(v, float) and np.isnan(v)):
                return None
            s = str(v).strip()
            if s.endswith(".0") and s[:-2].isdigit():
                return s[:-2]
            return s
        # add both nodes with attributes, then the edge
        for _, row in tqdm(
            primekg_df.iterrows(), total=len(primekg_df), desc="Building graph!"
        ):
            G.add_node(row["x_index"], name=row["x_name"], type=row["x_type"], kg_id=normalize_id(row["x_id"]))
            G.add_node(row["y_index"], name=row["y_name"], type=row["y_type"], kg_id=normalize_id(row["y_id"]))
            G.add_edge(row["x_index"], row["y_index"], relation=row["relation"])

        self.G = G
        # node sets by type, straight from the graph
        self.gene_set = {
            n for n, d in G.nodes(data=True) if d.get("type") == "gene/protein"
        }
        # diseases 
        diseases = {n for n, d in G.nodes(data=True) if d.get("type") == "disease"}
        # drugs
        drug = {n for n, d in G.nodes(data=True) if d.get("type") == "drug"}
        # effects
        effects = {
            n for n, d in G.nodes(data=True) if d.get("type") == "effect/phenotype"
        }
        # phenotype sets, all_phenotypes is their union
        self.phenotype_sets = {
            "gene/protein": self.gene_set,
            "diseases": diseases,
            "effects": effects,
            "drug":drug,
            "all_phenotypes": diseases
            | drug
            | effects
            | self.gene_set,
        }

        print(
            f"Graph built: {self.G.number_of_nodes()} nodes, {self.G.number_of_edges()} edges"
        )

    def save_cache(self, cache_file):
        """Pickle the graph and node sets to `cache_file`."""
        cache_data = {
            "graph": self.G,
            "gene_set": self.gene_set,
            "phenotype_sets": self.phenotype_sets,
        }
        with open(cache_file, "wb") as f:
            pickle.dump(cache_data, f)
        print(f"Cache saved to {cache_file}")

    def load_cache(self, cache_file):
        """Restore the graph and node sets from `cache_file`."""
        with open(cache_file, "rb") as f:
            cache_data = pickle.load(f)

        self.G = cache_data["graph"]
        self.gene_set = cache_data["gene_set"]
        self.phenotype_sets = cache_data["phenotype_sets"]
        print(f"Cache loaded from {cache_file}")

    # --> single metrics <-- #
    def calculate_shortest_paths(
        self, phenotypes, max_distance=4, exclude_node_types=None
    ):
        """Shortest path length from each phenotype to every gene within `max_distance`."""
        if self.G is None:
            raise ValueError("Graph not loaded")

        # drop the excluded node types on a copy
        if exclude_node_types is not None:
            nodes_to_exclude = {
                n
                for n, d in self.G.nodes(data=True)
                if d.get("type") in exclude_node_types
            }
            G_sub = self.G.copy()
            G_sub.remove_nodes_from(nodes_to_exclude)
        else:
            G_sub = self.G

        results = []

        for phenotype in tqdm(
            phenotypes, desc="Calculating distances"
        ):
            # BFS from the phenotype node, skip it when absent
            if phenotype not in G_sub:
                print(f"Phenotype {phenotype} not found in graph, skipping.")
                continue
            distances = nx.single_source_shortest_path_length(
                G_sub, phenotype, cutoff=max_distance
            )

            # keep gene nodes only
            for gene, distance in distances.items():
                if gene in self.gene_set and gene in G_sub:
                    results.append(
                        {
                            "phenotype_idx": phenotype,
                            "phenotype_name": G_sub.nodes[phenotype]["name"],
                            "gene_idx": gene,
                            "gene_name": G_sub.nodes[gene]["name"],
                            "graph_distance": distance,
                        }
                    )

        return pd.DataFrame(results)

    def calculate_adamic_adar(self, phenotypes, exclude_node_types=None):
        """Adamic-Adar index between each phenotype and every gene."""
        if self.G is None:
            raise ValueError("Graph not loaded")

        # drop the excluded node types on a copy
        if exclude_node_types is not None:
            nodes_to_exclude = {
                n
                for n, d in self.G.nodes(data=True)
                if d.get("type") in exclude_node_types
            }
            G_sub = self.G.copy()
            G_sub.remove_nodes_from(nodes_to_exclude)
            gene_set = {n for n in self.gene_set if n in G_sub}
        else:
            G_sub = self.G
            gene_set = self.gene_set

        results = []
        for phenotype in tqdm(phenotypes, desc="Calculating Adamic-Adar"):
            if phenotype not in G_sub:
                print(f"Phenotype {phenotype} not found in graph, skipping.")
                continue
            node_pairs = [(phenotype, gene) for gene in gene_set]

            # skip the phenotype when adamic_adar_index divides by zero
            try:
                aa_scores = nx.adamic_adar_index(G_sub, node_pairs)
                for phenotype_key, gene, score in aa_scores:
                    results.append(
                        {
                            "phenotype_idx": phenotype,
                            "phenotype_name": G_sub.nodes[phenotype]["name"],
                            "gene_idx": gene,
                            "gene_name": G_sub.nodes[gene]["name"],
                            "adamic_adar": score,
                        }
                    )
            except ZeroDivisionError:
                print(f"ZeroDivisionError encountered for phenotype {phenotype}, skipping.")
                continue
        return pd.DataFrame(results)

    def calculate_personalized_pagerank(
        self, phenotypes, alpha=0.85, exclude_node_types=None
    ):
        """Personalised PageRank from each phenotype, read off at every gene."""
        if self.G is None:
            raise ValueError("Graph not loaded")

        # drop the excluded node types on a copy
        if exclude_node_types is not None:
            nodes_to_exclude = {
                n
                for n, d in self.G.nodes(data=True)
                if d.get("type") in exclude_node_types
            }
            G_sub = self.G.copy()
            G_sub.remove_nodes_from(nodes_to_exclude)
            gene_set = {n for n in self.gene_set if n in G_sub}
        else:
            G_sub = self.G
            gene_set = self.gene_set

        results = []
        for phenotype in tqdm(phenotypes, desc="Calculating Personalized PageRank"):
            if phenotype not in G_sub:
                print(f"Phenotype {phenotype} not found in graph, skipping.")
                continue
            # restart at the phenotype node only
            personalization = {node: 0 for node in G_sub.nodes()}
            personalization[phenotype] = 1.0
            ppr_scores = nx.pagerank(
                G_sub, alpha=alpha, personalization=personalization
            )
            for gene in gene_set:
                score = ppr_scores.get(gene, 0)
                results.append(
                    {
                        "phenotype_idx": phenotype,
                        "phenotype_name": G_sub.nodes[phenotype]["name"],
                        "gene_idx": gene,
                        "gene_name": G_sub.nodes[gene]["name"],
                        "pagerank_score": score,
                    }
                )
        return pd.DataFrame(results)

    def score_phenotypes(self, phenotypes, max_distance=4, alpha=0.85):
        """Shortest path, Adamic-Adar and PageRank per (phenotype, gene), on the full graph."""
        if self.G is None:
            raise ValueError("Graph not loaded")

        # --> combined metrics <-- #
        distances_df = self.calculate_shortest_paths(phenotypes, max_distance)
        adamic_df = self.calculate_adamic_adar(phenotypes)
        pagerank_df = self.calculate_personalized_pagerank(phenotypes, alpha=alpha)

        # merge on phenotype and gene, unreached pairs get distance 5 and score 0
        results = distances_df.merge(
            adamic_df[
                [
                    "phenotype_idx",
                    "phenotype_name",
                    "gene_name",
                    "gene_idx",
                    "adamic_adar",
                ]
            ],
            on=["phenotype_idx", "phenotype_name", "gene_name", "gene_idx"],
            how="outer",
        ).fillna({"graph_distance": 5, "adamic_adar": 0})
        # add PageRank
        results = results.merge(
            pagerank_df[
                [
                    "phenotype_idx",
                    "phenotype_name",
                    "gene_name",
                    "gene_idx",
                    "pagerank_score",
                ]
            ],
            on=["phenotype_idx", "phenotype_name", "gene_name", "gene_idx"],
            how="outer",
        ).fillna({"graph_distance": 5, "pagerank_score": 0, "adamic_adar": 0})

        return pl.from_pandas(results)

    def score_phenotypes_with_and_without_exclusions(
        self, phenotypes, max_distance=4, alpha=0.85, exclude_node_types=None
    ):
        """Shortest path, Adamic-Adar and PageRank on the full graph and with `exclude_node_types` removed."""
        if self.G is None:
            raise ValueError("Graph not loaded")

        # --> full graph <-- #
        logger.info("Calculating Graph Scores with all nodes...")

        distances_df = self.calculate_shortest_paths(phenotypes, max_distance)
        adamic_df = self.calculate_adamic_adar(phenotypes)
        pagerank_df = self.calculate_personalized_pagerank(phenotypes, alpha=alpha)

        # merge on phenotype and gene, unreached pairs get distance 5 and score 0
        merge_keys = ["phenotype_idx", "phenotype_name", "gene_name", "gene_idx"]
        results = distances_df.merge(
            adamic_df[merge_keys + ["adamic_adar"]],
            on=merge_keys,
            how="outer",
        ).fillna({"graph_distance": 5, "adamic_adar": 0})
        results = results.merge(
            pagerank_df[merge_keys + ["pagerank_score"]],
            on=merge_keys,
            how="outer",
        ).fillna({"graph_distance": 5, "pagerank_score": 0, "adamic_adar": 0})

        # --> excluded node types removed <-- #
        logger.info("Calculating Graph Scores with GO Term exclusions...")
        distances_excl_df = self.calculate_shortest_paths(
            phenotypes, max_distance, exclude_node_types=exclude_node_types
        )
        adamic_excl_df = self.calculate_adamic_adar(
            phenotypes, exclude_node_types=exclude_node_types
        )
        pagerank_excl_df = self.calculate_personalized_pagerank(
            phenotypes, alpha=alpha, exclude_node_types=exclude_node_types
        )

        # suffix the score columns with _with_exclusions
        distances_excl_df = distances_excl_df.rename(
            columns={"graph_distance": "graph_distance_with_exclusions"}
        )
        adamic_excl_df = adamic_excl_df.rename(
            columns={"adamic_adar": "adamic_adar_with_exclusions"}
        )
        pagerank_excl_df = pagerank_excl_df.rename(
            columns={"pagerank_score": "pagerank_score_with_exclusions"}
        )

        # merge onto the full-graph scores
        results = results.merge(
            distances_excl_df[merge_keys + ["graph_distance_with_exclusions"]],
            on=merge_keys,
            how="outer",
        )
        results = results.merge(
            adamic_excl_df[merge_keys + ["adamic_adar_with_exclusions"]],
            on=merge_keys,
            how="outer",
        )
        results = results.merge(
            pagerank_excl_df[merge_keys + ["pagerank_score_with_exclusions"]],
            on=merge_keys,
            how="outer",
        )

        # same fill values as the full graph
        results = results.fillna(
            {
                "graph_distance_with_exclusions": 5,
                "adamic_adar_with_exclusions": 0,
                "pagerank_score_with_exclusions": 0,
            }
        )

        return pl.from_pandas(results)
    
    # --> TxGNN <-- #
    def attach_txgnn(self, export_dir, device="cpu"):
        """Attach the offline TxGNN scorer."""
        self.txgnn = TxGNNOfflineScorer(export_dir=export_dir)

    def txgnn_score_pair(self, idx_a, idx_b, rel=None):
        """DistMult score of two PrimeKG nodes, relation picked from their types when `rel` is None."""
        if not hasattr(self, "txgnn"):
                raise ValueError("Call attach_txgnn_numpy(export_dir) first")

        # node types and PrimeKG ids
        a = self.G.nodes[idx_a]
        b = self.G.nodes[idx_b]
        a_type, b_type = a.get("type"), b.get("type")
        a_id, b_id = a.get("kg_id"), b.get("kg_id")
        if a_id is None or b_id is None:
            raise ValueError("Missing kg_id on nodes; rebuild graph cache with kg_id for x and y.")
        # (source type, target type) -> TxGNN relation
        RELATION_MAP = {
        # self
        ("anatomy","anatomy"): "anatomy_anatomy",
        ("biological_process","biological_process"): "bioprocess_bioprocess",
        ("cellular_component","cellular_component"): "cellcomp_cellcomp",
        ("disease","disease"): "disease_disease",
        ("drug","drug"): "drug_drug",
        ("effect/phenotype","effect/phenotype"): "phenotype_phenotype",
        ("exposure","exposure"): "exposure_exposure",
        ("molecular_function","molecular_function"): "molfunc_molfunc",
        ("pathway","pathway"): "pathway_pathway",
        ("gene/protein","gene/protein"): "protein_protein",

        # disease <-> phenotype
        ("disease","effect/phenotype"): "disease_phenotype_positive",
        ("effect/phenotype","disease"): "rev_disease_phenotype_positive",

        # drug <-> disease
        ("drug","disease"): "indication",
        ("disease","drug"): "rev_indication",

        # drug <-> phenotype
        ("drug","effect/phenotype"): "drug_effect",
        ("effect/phenotype","drug"): "rev_drug_effect",

        # drug <-> gene
        ("drug","gene/protein"): "drug_protein",
        ("gene/protein","drug"): "rev_drug_protein",

        # disease <-> gene
        ("gene/protein","disease"): "disease_protein",
        ("disease","gene/protein"): "rev_disease_protein",

        # phenotype <-> gene
        ("gene/protein","effect/phenotype"): "phenotype_protein",
        ("effect/phenotype","gene/protein"): "rev_phenotype_protein",

        # anatomy <-> gene
        ("gene/protein","anatomy"): "anatomy_protein_present",
        ("anatomy","gene/protein"): "rev_anatomy_protein_present",

        # GO terms <-> gene
        ("gene/protein","biological_process"): "bioprocess_protein",
        ("biological_process","gene/protein"): "rev_bioprocess_protein",

        ("gene/protein","cellular_component"): "cellcomp_protein",
        ("cellular_component","gene/protein"): "rev_cellcomp_protein",

        ("gene/protein","molecular_function"): "molfunc_protein",
        ("molecular_function","gene/protein"): "rev_molfunc_protein",

        # pathway <-> gene
        ("gene/protein","pathway"): "pathway_protein",
        ("pathway","gene/protein"): "rev_pathway_protein",

        # exposure <-> others
        ("exposure","biological_process"): "exposure_bioprocess",
        ("biological_process","exposure"): "rev_exposure_bioprocess",

        ("exposure","cellular_component"): "exposure_cellcomp",
        ("cellular_component","exposure"): "rev_exposure_cellcomp",

        ("exposure","disease"): "exposure_disease",
        ("disease","exposure"): "rev_exposure_disease",

        ("exposure","molecular_function"): "exposure_molfunc",
        ("molecular_function","exposure"): "rev_exposure_molfunc",

        ("exposure","gene/protein"): "exposure_protein",
        ("gene/protein","exposure"): "rev_exposure_protein",
        }
        # None when the type pair has no relation
        if rel is None:
            rel = RELATION_MAP.get((a_type, b_type))
            if rel is None:
                return None

        return self.txgnn.distmult(a_type, a_id, rel, b_type, b_id)
                
    def add_txgnn_scores_to_pairs_df(
        self,
        pairs_df: pd.DataFrame,
        rel=None,
        out_col="txgnn_distmult",
    ):
        """Add `out_col` to a (phenotype_idx, gene_idx) frame, the DistMult score or NaN."""
        if not hasattr(self, "txgnn"):
            raise ValueError("Call attach_txgnn(export_dir) first")

        # start at NaN
        scores = np.empty(len(pairs_df), dtype=np.float32)
        scores[:] = np.nan

        it = zip(
            pairs_df["phenotype_idx"].to_numpy(),
            pairs_df["gene_idx"].to_numpy(),
        )

        # score each pair, any failure leaves NaN
        for i, (p, g) in enumerate(
            tqdm(it, total=len(pairs_df), desc=f"TxGNN -> {out_col}")
        ):
            try:
                s = self.txgnn_score_pair(int(p), int(g), rel=rel)
                if s is not None:
                    scores[i] = float(s)
            except Exception:
                pass

        pairs_df[out_col] = scores
        return pairs_df


################################################
####                  setup                 ####
################################################


def setup_graph_scorer(primekg_path, cache_file_path):
    """GraphScorer loaded from the cache, or built from PrimeKG `kg.csv` and cached."""
    graph_scorer = GraphScorer()

    # load the cache, build and save it otherwise
    if os.path.exists(cache_file_path):
        graph_scorer.load_cache(cache_file_path)
    else:
        print(f"Building the PrimeKG graph, first run only, cached to {cache_file_path} afterwards")
        primekg_df = pd.read_csv(f"{primekg_path}/kg.csv", low_memory=False)
        graph_scorer.build_graph(primekg_df)
        graph_scorer.save_cache(cache_file_path)

    return graph_scorer


def score_phenotypes(graph_scorer, phenotypes, max_distance=4, alpha=0.85):
    """Score phenotypes with a GraphScorer."""
    if not isinstance(graph_scorer, GraphScorer):
        raise ValueError("graph_scorer must be an instance of GraphScorer")

    return graph_scorer.analyze_phenotypes(phenotypes, max_distance)


################################################
####          variant table scores          ####
################################################


def apply_graph_scores(variant_table: pl.DataFrame, config):
    """Attach graph-based scores to the variant table and return the scorer."""
    # --> scorer and maps <-- #
    graph_scorer = setup_graph_scorer(
        config.graph_scorer.primekg_dir, config.graph_scorer.cache_file
    )
    # attach the offline TxGNN export
    graph_scorer.attach_txgnn(
        export_dir=config.graph_scorer.txgnn_export_dir,
        device="cpu",
    )
    # locus to phenotype node, Ensembl id to PrimeKG gene node
    marker_id_phenotype_map_df = pl.read_csv(
        config.dataset.gwas_endpoint_file, separator="\t"
    )
    # a repeated row would repeat every gene row of its locus in the join below
    marker_id_phenotype_map_df = marker_id_phenotype_map_df.unique(
        subset=["locus_id", "phenotype_idx"], maintain_order=True
    )
    ensembl_gene_ids_to_primekg_idx_map = pl.read_csv(
        config.graph_scorer.ensembl_to_primekg_idx_file, separator="\t"
    ).rename({"node_index": "gene_idx"})

    # --> phenotype nodes as int <-- #
    list_of_phenotypes = (
        marker_id_phenotype_map_df["phenotype_idx"].unique().to_list()
    )

    def safe_int(x):
        """int(x), None when it does not parse."""
        try:
            return int(x)
        except (ValueError, TypeError):
            return None

    list_of_phenotypes_int = [safe_int(x) for x in list_of_phenotypes]
    list_of_phenotypes = [x for x in list_of_phenotypes_int if x is not None]
    if len(list_of_phenotypes) != len(list_of_phenotypes_int):
        logger.warning(
            "Some phenotype indices could not be converted to int and were skipped."
        )

    # --> graph scores <-- #
    graph_scoring_results_df = graph_scorer.score_phenotypes_with_and_without_exclusions(
        list_of_phenotypes,
        max_distance=config.graph_scorer.max_distance,
        alpha=config.graph_scorer.alpha,
        exclude_node_types=config.graph_scorer.exclude_node_types,
    )

    # --> TxGNN scores per (phenotype, gene) <-- #
    tmp = graph_scoring_results_df.to_pandas()

    # node ids must be int for the graph lookup
    tmp["phenotype_idx"] = tmp["phenotype_idx"].astype(int)
    tmp["gene_idx"] = tmp["gene_idx"].astype(int)

    # rel=None routes each pair to the relation of its node types
    tmp = graph_scorer.add_txgnn_scores_to_pairs_df(
        tmp, rel=None, out_col="txgnn_distmult"
    )

    graph_scoring_results_df = pl.from_pandas(tmp)

    # --> map to locus and Ensembl gene <-- #
    graph_scoring_results_df = graph_scoring_results_df.with_columns(
        pl.col("phenotype_idx").cast(pl.String)
    )
    marker_id_phenotype_map_df = marker_id_phenotype_map_df.with_columns(
        pl.col("phenotype_idx").cast(pl.String)
    )
    graph_scoring_results_df = graph_scoring_results_df.join(
        marker_id_phenotype_map_df[["locus_id", "phenotype_idx"]],
        on="phenotype_idx",
        how="left",
    )
    graph_scoring_results_df = graph_scoring_results_df.join(
        ensembl_gene_ids_to_primekg_idx_map[["gene_idx", "ensembl_gene_id"]],
        on="gene_idx",
        how="left",
    )
    graph_scoring_results_df = graph_scoring_results_df.drop_nulls()

    # --> join onto the variant table <-- #
    variant_table = variant_table.join(
        graph_scoring_results_df[
            [
                "locus_id",
                "ensembl_gene_id",
                "graph_distance",
                "adamic_adar",
                "pagerank_score",
                "graph_distance_with_exclusions",
                "adamic_adar_with_exclusions",
                "pagerank_score_with_exclusions",
                "txgnn_distmult",
            ]
        ],
        left_on=["locus_id", "ensembl_gene_id"],
        right_on=["locus_id", "ensembl_gene_id"],
        how="left",
    )
    # unscored pairs get distance 5 and score 0
    variant_table = variant_table.with_columns(
        pl.col("graph_distance").fill_null(5),
        pl.col("adamic_adar").fill_null(0),
        pl.col("pagerank_score").fill_null(0),
        pl.col("graph_distance_with_exclusions").fill_null(5),
        pl.col("adamic_adar_with_exclusions").fill_null(0),
        pl.col("pagerank_score_with_exclusions").fill_null(0),
    )

    return variant_table, graph_scorer
