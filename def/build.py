import json
from pathlib import Path
from collections import defaultdict, deque
import pandas as pd


# CONFIG PERCORSI 

HERE = Path(__file__).resolve().parent

BASE = HERE.parent / "csv"

GO_OBO = BASE / "go-basic.obo"
HPO_OBO = BASE / "hp.obo"

GENE_HPO = BASE / "gene_hpo_matrix_binary_directOnly_namespace_Phenotypic abnormality.csv"


# PARSER GO OBO 
def parse_go_obo(obo_path: Path):
    terms = {}                      # estrae id, name, namespace
    parents = defaultdict(set)      # estrae is_a e part_of
    rel_partof = set()              # memorizza relazioni part_of

    with open(obo_path, encoding="utf-8") as f:
        current = None

        # parsing riga per riga
        for line in f:
            line = line.strip() # rimuovi spazi

            # quando trova [Term], inizia nuovo termine
            if line == "[Term]":
                current = {}

            # quando trova riga vuota, termina termine corrente
            elif line == "" and current:
                # se ha id, salva termine
                if "id" in current:
                    terms[current["id"]] = current
                # poi resetta current
                current = None

            # se sta leggendo un termine, parsa le info
            elif current is not None:
                if line.startswith("id: "):
                    current["id"] = line.split("id: ")[1]

                elif line.startswith("name: "):
                    current["name"] = line.split("name: ")[1]

                elif line.startswith("namespace: "):
                    current["namespace"] = line.split("namespace: ")[1]

                # estrae il parent per is_a
                elif line.startswith("is_a: "):
                    parent = line.split("is_a: ")[1].split(" ! ")[0]
                    parents[current.get("id", "")].add(parent)

                # estrae il parent per part_of
                elif line.startswith("relationship: part_of "):
                    parent = line.split("relationship: part_of ")[1].split(" ! ")[0]
                    parents[current.get("id", "")].add(parent)
                    rel_partof.add((current.get("id", ""), parent))

        # flush ultimo term se file non finisce con riga vuota
        if current and "id" in current:
            terms[current["id"]] = current

    return terms, parents, rel_partof

# data una radice e una mappa child->parents, calcola depth
def compute_depths_from_root(root_id: str, parents_map: dict):
    children = defaultdict(list)
    for child, ps in parents_map.items():
        for p in ps:
            children[p].append(child)

    depth = {}
    depth[root_id] = 0
    q = deque([root_id])

    # BFS per calcolare depth
    while q:
        node = q.popleft()
        for ch in children.get(node, []):
            if ch not in depth:
                depth[ch] = depth[node] + 1
                q.append(ch)

    return depth



# PARSER HPO OBO
def parse_hpo_obo(obo_path: Path):
    terms = {}
    parents = defaultdict(set)

    # parsing riga per riga
    with open(obo_path, encoding="utf-8") as f:
        current = None

        for line in f:
            line = line.strip()

            if line == "[Term]":
                if current and "id" in current:
                    terms[current["id"]] = current
                current = {}

            elif current is not None:
                if line.startswith("id: "):
                    current["id"] = line.split("id: ")[1]
                elif line.startswith("name: "):
                    current["name"] = line.split("name: ")[1]
                elif line.startswith("is_a: "):
                    parent = line.split("is_a: ")[1].split(" ! ")[0]
                    parents[current["id"]].add(parent)

        if current and "id" in current:
            terms[current["id"]] = current

    # profondità HPO
    ROOT = "HP:0000001"
    depth = compute_depths_from_root(ROOT, parents)

    return terms, parents, depth


# CARICAMENTI MATRICI GENE-GO
def load_gene_go_matrix(csv_path: Path):
    df = pd.read_csv(csv_path, sep=";", low_memory=False)

    first = df.columns[0]
    if first != "ENTREZID":
        df.rename(columns={first: "ENTREZID"}, inplace=True)

    # riga "go_depth" va rimossa (dati di profondità non presenti per tutti i go)
    df = df[df["ENTREZID"].astype(str).str.lower() != "go_depth"]

    genes = df["ENTREZID"].dropna().astype(str).unique().tolist()
    df = df[df["ENTREZID"].astype(str).isin(genes)]

    # normalizza GO. -> GO: per matchare con obo
    go_cols = [c for c in df.columns if c.startswith("GO.")]
    mapping = {c: c.replace("GO.", "GO:") for c in go_cols}
    df.rename(columns=mapping, inplace=True)

    return df, list(mapping.values()), genes


# CARICAMENTO MATRICE GENE-HPO
def load_gene_hpo(genes, gene_hpo_path: Path):
    df = pd.read_csv(gene_hpo_path, sep=";", low_memory=False)

    first = df.columns[0]
    if first != "ENTREZID":
        df.rename(columns={first: "ENTREZID"}, inplace=True)

    df = df[df["ENTREZID"].astype(str).isin(genes)]

    # normalizza HP. -> HP: per matchare con obo
    hpo_cols_raw = [c for c in df.columns if c.startswith("HP:") or c.startswith("HP.")]
    hpo_mapping = {c: c.replace("HP.", "HP:") for c in hpo_cols_raw}
    df.rename(columns=hpo_mapping, inplace=True)

    cols_hpo = list(hpo_mapping.values())
    return df, cols_hpo


# GENERAZIONE SINGOLO JSON (BP / CC / MF)
def build_graph_for_namespace(
    *,
    go_terms,
    go_parents,
    go_rel_partof,
    hpo_terms,
    hpo_parents,
    hpo_depth,
    go_root: str,
    target_namespace: str,
    gene_go_csv: Path,
    output_json: Path,
):
    # depth GO specifica per root (BP/CC/MF)
    go_depth = compute_depths_from_root(go_root, go_parents)

    # carica matrici gene-go e gene-hpo (filtrata sui geni presenti nel gene-go)
    df_go, cols_go, genes = load_gene_go_matrix(gene_go_csv)
    df_hpo, cols_hpo = load_gene_hpo(genes, GENE_HPO)

    # NODI
    gene_nodes = [{"id": g, "type": "gene"} for g in genes]

    # Prende i GO dell'ontologia GO OBO che sono nel namespace target
    go_nodes = []
    for go_id, info in go_terms.items():
        if info.get("namespace") == target_namespace:
            go_nodes.append(
                {
                    "id": go_id,
                    "name": info.get("name", ""),
                    "type": "go",
                    "namespace": target_namespace,
                    "depth": go_depth.get(go_id),
                }
            )
    go_set = {n["id"] for n in go_nodes}

    # Prende gli HPO (tutti)
    hpo_nodes = []
    for hpo_id, info in hpo_terms.items():
        hpo_nodes.append(
            {
                "id": hpo_id,
                "name": info.get("name", ""),
                "type": "hpo",
                "depth": hpo_depth.get(hpo_id),
            }
        )

    # ARCHI
    # GO -> GO (solo dentro namespace target)
    edges_go = []
    for child, ps in go_parents.items():
        if child not in go_set:
            continue
        for parent in ps:
            if parent not in go_set:
                continue
            rel = "part_of" if (child, parent) in go_rel_partof else "is_a"
            edges_go.append({"source": child, "target": parent, "rel": rel})

    # HPO -> HPO (tutti)
    edges_hpo = []
    for child, ps in hpo_parents.items():
        for parent in ps:
            edges_hpo.append({"source": child, "target": parent, "rel": "is_a"})

    # gene -> GO
    edges_gene_go = []
    for _, row in df_go.iterrows():
        gid = str(row["ENTREZID"])
        for go_id in cols_go:
            try:
                val = float(row[go_id])
            except Exception:
                continue
            if val > 0 and go_id in go_set:
                edges_gene_go.append({"source": gid, "target": go_id, "rel": "gene_go"})

    # gene -> HPO
    edges_gene_hpo = []
    for _, row in df_hpo.iterrows():
        gid = str(row["ENTREZID"])
        for hpo_id in cols_hpo:
            try:
                val = float(row[hpo_id])
            except Exception:
                continue
            if val > 0:
                edges_gene_hpo.append({"source": gid, "target": hpo_id, "rel": "gene_hpo"})

    # scrive JSON finale, con due liste: nodi ed archi
    graph = {
        "nodes": gene_nodes + go_nodes + hpo_nodes,
        "edges": edges_go + edges_hpo + edges_gene_go + edges_gene_hpo,
    }

    output_json.write_text(json.dumps(graph, indent=2), encoding="utf-8")

    # log di riepilogo
    print(f"✅ GENERATO: {output_json}")
    print(f"→ GENI: {len(gene_nodes)}")
    print(f"→ GO ({target_namespace}): {len(go_nodes)}")
    print(f"→ HPO: {len(hpo_nodes)}")
    print(f"→ edges GO: {len(edges_go)}")
    print(f"→ edges HPO: {len(edges_hpo)}")
    print(f"→ edges gene→GO: {len(edges_gene_go)}")
    print(f"→ edges gene→HPO: {len(edges_gene_hpo)}")
    print("-" * 60)


def main():
    # parse una volta sola
    go_terms, go_parents, go_rel_partof = parse_go_obo(GO_OBO)
    hpo_terms, hpo_parents, hpo_depth = parse_hpo_obo(HPO_OBO)

    # configurazione dei 3 casi
    jobs = [
        {
            "label": "BP",
            "go_root": "GO:0008150",
            "target_namespace": "biological_process",
            "gene_go_csv": BASE / "gene_go_matrix_propF_rel-is_a-part_of_ont-BP_withDepths.csv",
            "output_json": Path("graph_go_bp_gene_hpo.json"),
        },
        {
            "label": "CC",
            "go_root": "GO:0005575",
            "target_namespace": "cellular_component",
            "gene_go_csv": BASE / "gene_go_matrix_propF_rel-is_a-part_of_ont-CC_withDepths.csv",
            "output_json": Path("graph_go_cc_gene_hpo.json"),
        },
        {
            "label": "MF",
            "go_root": "GO:0003674",
            "target_namespace": "molecular_function",
            "gene_go_csv": BASE / "gene_go_matrix_propF_rel-is_a-part_of_ont-MF_withDepths.csv",
            "output_json": Path("graph_go_mf_gene_hpo.json"),
        },
    ]

    for job in jobs:
        print(f"=== BUILD {job['label']} ===")
        build_graph_for_namespace(
            go_terms=go_terms,
            go_parents=go_parents,
            go_rel_partof=go_rel_partof,
            hpo_terms=hpo_terms,
            hpo_parents=hpo_parents,
            hpo_depth=hpo_depth,
            go_root=job["go_root"],
            target_namespace=job["target_namespace"],
            gene_go_csv=job["gene_go_csv"],
            output_json=job["output_json"],
        )


if __name__ == "__main__":
    main()
