import glob
import os
import re
import shutil
import threading
import time
import uuid
from collections import defaultdict

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.nn import GATConv, GCNConv, global_mean_pool

# Custom preprocessing used by the existing website.
from FilterAtomsGraphs import create_graphs
from vina_docking import DockingConfig, DockingError, dock_cholesterol


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

ATOM_TYPE_NAMES = (
    "C", "CA", "CB", "CD", "CD1", "CD2", "CE", "CE1", "CE2", "CE3",
    "CG", "CG1", "CG2", "CH2", "CZ", "CZ2", "CZ3", "O", "OH", "OD1",
    "OD2", "OE1", "OE2", "OG", "OG1", "N", "NE", "NE1", "NE2", "ND1",
    "ND2", "NZ", "NH1", "NH2", "SD", "SG", "UNKNOWN",
)

PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "ASX", "GLX", "SEC", "PYL", "MSE",
}

_CLR_TAG_RE = re.compile(r"_CLR(-?\d+)([A-Za-z0-9])(?=[_-])")


def parse_clr_from_path(path):
    """Extract ``(residue_number, chain_id)`` from a generated site filename."""
    match = _CLR_TAG_RE.search(os.path.basename(path))
    if not match:
        return None
    return int(match.group(1)), match.group(2)


def robust_rmtree(path, retries=5, delay=0.2):
    if not os.path.exists(path):
        return
    for _ in range(retries):
        try:
            shutil.rmtree(path)
            return
        except OSError:
            time.sleep(delay)
    print(f"Warning: Could not fully delete temp dir {path} due to file locks.")


class GAT(torch.nn.Module):
    def __init__(self, in_channels, out_channels, dropout_p=0.1):
        super().__init__()
        self.gat = GATConv(
            in_channels, out_channels, heads=1, concat=True, edge_dim=1
        )
        self.pool = global_mean_pool
        self.dropout = nn.Dropout(p=dropout_p)
        self.norm = nn.BatchNorm1d(out_channels)
        self.linear = torch.nn.Linear(out_channels, 1)

    def forward(self, x, edge_index, edge_attr, batch):
        out, attn_weights = self.gat(
            x,
            edge_index,
            edge_attr,
            return_attention_weights=True,
        )
        out = self.dropout(out)
        out = self.pool(out, batch)
        out = self.norm(out)
        out = self.dropout(out)
        out = self.linear(out)
        return out, attn_weights


class GCN(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.conv1 = GCNConv(input_dim, 32)
        self.bn1 = nn.BatchNorm1d(32)
        self.conv2 = GCNConv(32, 64)
        self.bn2 = nn.BatchNorm1d(64)
        self.conv3 = GCNConv(64, 128)
        self.bn3 = nn.BatchNorm1d(128)
        self.dropout_gcn = nn.Dropout(0.2)
        self.dropout = nn.Dropout(0.6)
        self.fc1 = nn.Linear(128, 64)
        self.out = nn.Linear(64, 1)

    def forward(self, data):
        x = data.x
        edge_index = data.edge_index
        edge_weight = data.edge_attr
        batch = data.batch

        x = self.dropout_gcn(F.relu(self.bn1(self.conv1(x, edge_index, edge_weight))))
        x = self.dropout_gcn(F.relu(self.bn2(self.conv2(x, edge_index, edge_weight))))
        x = self.dropout_gcn(F.relu(self.bn3(self.conv3(x, edge_index, edge_weight))))
        x = global_mean_pool(x, batch)
        x = self.dropout(F.relu(self.fc1(x)))
        return self.out(x)


class CNN2D(nn.Module):
    def __init__(self, input_channels):
        super().__init__()
        self.conv1 = nn.Sequential(
            nn.Conv2d(input_channels, 32, kernel_size=3, padding=1),
            nn.BatchNorm2d(32),
            nn.ReLU(),
        )
        self.pool1 = nn.MaxPool2d(2, 2)
        self.conv2 = nn.Sequential(
            nn.Conv2d(32, 64, kernel_size=3, padding=1),
            nn.BatchNorm2d(64),
            nn.ReLU(),
        )
        self.pool2 = nn.MaxPool2d(2, 2)
        self.conv3 = nn.Sequential(
            nn.Conv2d(64, 128, kernel_size=3, padding=1),
            nn.BatchNorm2d(128),
            nn.ReLU(),
        )
        self.pool3 = nn.MaxPool2d(2, 2)
        self.flatten = nn.Flatten()
        self.fc1 = nn.Linear(128 * 4 * 18, 128)
        self.dropout = nn.Dropout(0.5)
        self.fc2 = nn.Linear(128, 1)

    def forward(self, x):
        x = self.pool1(self.conv1(x))
        x = self.pool2(self.conv2(x))
        x = self.pool3(self.conv3(x))
        x = self.flatten(x)
        x = self.dropout(torch.relu(self.fc1(x)))
        return self.fc2(x)


def _load_graph_payload(file_path):
    payload = np.load(file_path, allow_pickle=True).item()
    if not {"inverse_distance", "encoded_matrix"}.issubset(payload):
        raise ValueError(f"Graph file is missing required arrays: {file_path}")
    return payload


def organize_graph_gat(file_path, label=0):
    payload = _load_graph_payload(file_path)
    x = torch.tensor(payload["encoded_matrix"], dtype=torch.float32)
    adj = torch.tensor(payload["inverse_distance"], dtype=torch.float32)
    adj = adj / (adj.sum(dim=1, keepdim=True) + 1e-8)
    edge_index = (adj > 0).nonzero(as_tuple=False).t()
    edge_weight = adj[adj > 0]
    y = torch.tensor([label], dtype=torch.float32)
    return Data(x=x, edge_index=edge_index, edge_attr=edge_weight, y=y)


def organize_graph_gcn(file_path, label=0):
    payload = _load_graph_payload(file_path)
    x = torch.tensor(payload["encoded_matrix"], dtype=torch.float32)
    adj = torch.tensor(payload["inverse_distance"], dtype=torch.float32)
    edge_index = (adj > 0).nonzero(as_tuple=False).t()
    edge_weight = adj[adj > 0]
    y = torch.tensor([label], dtype=torch.float32)
    return Data(x=x, edge_index=edge_index, edge_attr=edge_weight, y=y)


def _decode_atom_type(one_hot_row):
    row = np.asarray(one_hot_row)
    if row.size == 0 or not np.any(row):
        return "UNKNOWN"
    index = int(np.argmax(row))
    return ATOM_TYPE_NAMES[index] if index < len(ATOM_TYPE_NAMES) else "UNKNOWN"


def _parse_pdb_atoms(pdb_path):
    """Parse first-model PDB atoms in file order, preserving PDB serial IDs."""
    atoms = []
    saw_model = False
    in_first_model = True

    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip().upper()
            if record == "MODEL":
                if saw_model:
                    break
                saw_model = True
                in_first_model = True
                continue
            if record == "ENDMDL" and saw_model:
                break
            if not in_first_model or record not in {"ATOM", "HETATM"}:
                continue

            try:
                atom_name = line[12:16].strip()
                element = line[76:78].strip().upper()
                if not element:
                    letters = "".join(ch for ch in atom_name if ch.isalpha())
                    element = letters[:1].upper()

                atoms.append(
                    {
                        "record": record,
                        "serial": int(line[6:11].strip()),
                        "atom_name": atom_name,
                        "alt_loc": line[16:17].strip(),
                        "residue_name": line[17:20].strip().upper(),
                        "chain_id": line[21:22].strip(),
                        "residue_number": int(line[22:26].strip()),
                        "insertion_code": line[26:27].strip(),
                        "x": float(line[30:38]),
                        "y": float(line[38:46]),
                        "z": float(line[46:54]),
                        "element": element,
                    }
                )
            except (TypeError, ValueError):
                continue

    return atoms


def _residue_key(atom):
    return (
        atom["chain_id"],
        atom["residue_number"],
        atom["insertion_code"],
        atom["residue_name"],
    )


def _squared_distance(atom_a, atom_b):
    return (
        (atom_a["x"] - atom_b["x"]) ** 2
        + (atom_a["y"] - atom_b["y"]) ** 2
        + (atom_a["z"] - atom_b["z"]) ** 2
    )


def _atom_type_match_fraction(atoms, encoded_matrix):
    comparable = 0
    matches = 0
    for atom, row in zip(atoms, encoded_matrix):
        expected = _decode_atom_type(row)
        if expected == "UNKNOWN":
            continue
        comparable += 1
        if atom["atom_name"].upper() == expected:
            matches += 1
    return matches / comparable if comparable else 1.0


def _site_atom_candidates(pdb_path, clr_key, cutoff=5.0):
    atoms = _parse_pdb_atoms(pdb_path)
    protein = [
        atom
        for atom in atoms
        if atom["element"] != "H"
        and (atom["record"] == "ATOM" or atom["residue_name"] in PROTEIN_RESIDUES)
    ]

    residue_number, chain_id = clr_key
    ligand = [
        atom
        for atom in atoms
        if atom["residue_number"] == residue_number
        and atom["chain_id"] == chain_id
        and atom["residue_name"] in {"CLR", "CHL"}
        and atom["element"] != "H"
    ]
    if not ligand:
        ligand = [
            atom
            for atom in atoms
            if atom["residue_number"] == residue_number
            and atom["chain_id"] == chain_id
            and atom["record"] == "HETATM"
            and atom["element"] != "H"
        ]

    candidates = [("all protein atoms", protein)]
    if ligand:
        cutoff_sq = cutoff * cutoff
        atomwise = [
            atom
            for atom in protein
            if any(_squared_distance(atom, lig_atom) <= cutoff_sq for lig_atom in ligand)
        ]
        near_residue_keys = {_residue_key(atom) for atom in atomwise}
        residuewise = [
            atom for atom in protein if _residue_key(atom) in near_residue_keys
        ]
        candidates.insert(0, (f"protein atoms within {cutoff:g} A", atomwise))
        candidates.insert(1, (f"protein residues within {cutoff:g} A", residuewise))

    return candidates


def _map_graph_nodes_to_pdb(
    original_pdb_path,
    clr_key,
    encoded_matrix,
    generated_pdb_paths=None,
):
    """
    Match graph nodes to PDB records only when atom count and atom-name order agree.

    Returning ``None`` is intentional: it prevents the UI from displaying a
    plausible-looking but incorrect chain/residue/serial identifier.
    """
    expected_count = int(encoded_matrix.shape[0])
    candidates = []

    for label, atoms in _site_atom_candidates(original_pdb_path, clr_key):
        candidates.append((label, atoms))

    for generated_path in generated_pdb_paths or []:
        tagged_clr = parse_clr_from_path(generated_path)
        if tagged_clr is not None and tagged_clr != clr_key:
            continue
        parsed = _parse_pdb_atoms(generated_path)
        generated_protein = [
            atom
            for atom in parsed
            if atom["element"] != "H"
            and (atom["record"] == "ATOM" or atom["residue_name"] in PROTEIN_RESIDUES)
        ]
        candidates.insert(0, ("generated filtered PDB", generated_protein))

    best = None
    seen_serial_orders = set()
    for label, atoms in candidates:
        serial_order = tuple(atom["serial"] for atom in atoms)
        if serial_order in seen_serial_orders:
            continue
        seen_serial_orders.add(serial_order)
        if len(atoms) != expected_count:
            continue
        match_fraction = _atom_type_match_fraction(atoms, encoded_matrix)
        if best is None or match_fraction > best[0]:
            best = (match_fraction, label, atoms)

    if best is None:
        return None, (
            f"Could not safely align {expected_count} graph nodes with PDB atoms; "
            "node-level importance is available without PDB identifiers."
        )
    if best[0] < 0.85:
        return None, (
            "Graph/PDB atom order validation was below 85%; identifiers were hidden "
            "to avoid highlighting the wrong atoms."
        )

    return best[2], None


def _positive_normalize(values):
    values = np.asarray(values, dtype=np.float64)
    values = np.where(np.isfinite(values), values, 0.0)
    values = np.maximum(values, 0.0)
    maximum = float(values.max()) if values.size else 0.0
    if maximum <= 1e-12:
        return np.zeros_like(values, dtype=np.float64)
    return values / maximum


def _summarize_categories(records, key_name, output_name):
    grouped = defaultdict(lambda: {"count": 0, "importance_sum": 0.0})
    for record in records:
        label = record.get(key_name)
        if not label:
            continue
        grouped[label]["count"] += 1
        grouped[label]["importance_sum"] += float(record.get("importance", 0.0))

    summary = []
    for label, values in grouped.items():
        count = values["count"]
        summary.append(
            {
                output_name: label,
                "count": count,
                "mean_importance": values["importance_sum"] / max(count, 1),
            }
        )
    summary.sort(key=lambda item: (-item["count"], -item["mean_importance"], item[output_name]))
    return summary


def _build_atom_interpretation(
    node_scores,
    encoded_matrix,
    atom_metadata,
    method,
    include_identifiers=True,
    top_k=5,
    warning=None,
    top_edges=None,
):
    node_scores = np.asarray(node_scores, dtype=np.float64)[: encoded_matrix.shape[0]]
    normalized = _positive_normalize(node_scores)
    ranked_indices = [
        int(index)
        for index in np.argsort(normalized)[::-1]
        if normalized[index] > 0
    ][:top_k]

    detailed_atoms = []
    for rank, node_index in enumerate(ranked_indices, start=1):
        item = {
            "rank": rank,
            "node_index": node_index,
            "atom_type": _decode_atom_type(encoded_matrix[node_index]),
            "importance": float(normalized[node_index]),
        }
        if atom_metadata is not None and node_index < len(atom_metadata):
            atom = atom_metadata[node_index]
            item.update(
                {
                    "atom_serial": atom["serial"],
                    "atom_name": atom["atom_name"],
                    "residue_name": atom["residue_name"],
                    "residue_number": atom["residue_number"],
                    "chain_id": atom["chain_id"],
                    "insertion_code": atom["insertion_code"],
                }
            )
        detailed_atoms.append(item)

    detailed_residues = []
    if atom_metadata is not None:
        residue_scores = defaultdict(lambda: {"score": 0.0, "atom_count": 0, "atom": None})
        for node_index, score in enumerate(normalized):
            if score <= 0 or node_index >= len(atom_metadata):
                continue
            atom = atom_metadata[node_index]
            key = _residue_key(atom)
            residue_scores[key]["score"] += float(score)
            residue_scores[key]["atom_count"] += 1
            residue_scores[key]["atom"] = atom

        ordered = sorted(
            residue_scores.values(),
            key=lambda item: -item["score"],
        )[:top_k]
        maximum = max((item["score"] for item in ordered), default=0.0)
        for rank, values in enumerate(ordered, start=1):
            atom = values["atom"]
            detailed_residues.append(
                {
                    "rank": rank,
                    "residue_name": atom["residue_name"],
                    "residue_number": atom["residue_number"],
                    "chain_id": atom["chain_id"],
                    "insertion_code": atom["insertion_code"],
                    "importance": values["score"] / maximum if maximum > 0 else 0.0,
                    "important_atom_count": values["atom_count"],
                }
            )

    atom_type_summary = _summarize_categories(detailed_atoms, "atom_type", "atom_type")
    residue_type_summary = _summarize_categories(
        detailed_atoms,
        "residue_name",
        "residue_name",
    )

    if not ranked_indices and warning is None:
        warning = "No positive importance signal was produced for this site."
    if atom_metadata is None and warning is None:
        warning = "PDB identifiers are unavailable because graph/PDB atom alignment could not be validated."

    result = {
        "method": method,
        "mapping_status": "matched" if atom_metadata is not None else "unavailable",
        "top_atom_types": atom_type_summary,
        "top_residue_types": residue_type_summary,
        "warning": warning,
    }
    if include_identifiers:
        result["top_atoms"] = detailed_atoms
        result["top_residues"] = detailed_residues
        if top_edges is not None:
            result["top_edges"] = top_edges
    return result


def _edge_based_interpretation(
    edge_runs,
    encoded_matrix,
    atom_metadata,
    method,
    include_identifiers,
    mapping_warning,
    top_k=5,
):
    """Aggregate signed positive edge evidence across all ensemble models."""
    across_runs = defaultdict(list)
    for edge_index, importance in edge_runs:
        per_run = defaultdict(list)
        sources = edge_index[0].reshape(-1)
        targets = edge_index[1].reshape(-1)
        values = np.asarray(importance).reshape(-1)
        for source, target, value in zip(sources, targets, values):
            source = int(source)
            target = int(target)
            if source == target or source >= encoded_matrix.shape[0] or target >= encoded_matrix.shape[0]:
                continue
            value = float(value)
            if not np.isfinite(value) or value <= 0:
                continue
            key = tuple(sorted((source, target)))
            per_run[key].append(value)
        for key, values_for_edge in per_run.items():
            across_runs[key].append(float(np.mean(values_for_edge)))

    candidates = []
    for (source, target), values in across_runs.items():
        raw_importance = float(np.mean(values))
        candidate = {
            "source": source,
            "target": target,
            "raw_importance": raw_importance,
            "model_frequency": len(values),
        }
        if atom_metadata is not None:
            source_atom = atom_metadata[source]
            target_atom = atom_metadata[target]
            if _residue_key(source_atom) == _residue_key(target_atom):
                continue
            distance = _squared_distance(source_atom, target_atom) ** 0.5
            if distance < 4.0:
                continue
            candidate["distance_angstrom"] = distance
            candidate["residue_pair"] = tuple(
                sorted((_residue_key(source_atom), _residue_key(target_atom)))
            )
        candidates.append(candidate)

    # Preserve the notebook's idea of one representative interaction per
    # residue pair, but use exact residue IDs for a single uploaded structure.
    if atom_metadata is not None:
        best_by_pair = {}
        for candidate in candidates:
            pair = candidate["residue_pair"]
            current = best_by_pair.get(pair)
            if current is None or candidate["raw_importance"] > current["raw_importance"]:
                best_by_pair[pair] = candidate
        candidates = list(best_by_pair.values())

    candidates.sort(key=lambda item: -item["raw_importance"])
    chosen = candidates[:top_k]
    maximum = max((item["raw_importance"] for item in chosen), default=0.0)

    node_scores = np.zeros(encoded_matrix.shape[0], dtype=np.float64)
    top_edges = []
    for rank, candidate in enumerate(chosen, start=1):
        normalized = candidate["raw_importance"] / maximum if maximum > 0 else 0.0
        source = candidate["source"]
        target = candidate["target"]
        node_scores[source] += normalized
        node_scores[target] += normalized

        if include_identifiers and atom_metadata is not None:
            source_atom = atom_metadata[source]
            target_atom = atom_metadata[target]
            top_edges.append(
                {
                    "rank": rank,
                    "importance": normalized,
                    "distance_angstrom": candidate["distance_angstrom"],
                    "model_frequency": candidate["model_frequency"],
                    "source": {
                        "node_index": source,
                        "atom_serial": source_atom["serial"],
                        "atom_name": source_atom["atom_name"],
                        "atom_type": _decode_atom_type(encoded_matrix[source]),
                        "residue_name": source_atom["residue_name"],
                        "residue_number": source_atom["residue_number"],
                        "chain_id": source_atom["chain_id"],
                        "insertion_code": source_atom["insertion_code"],
                    },
                    "target": {
                        "node_index": target,
                        "atom_serial": target_atom["serial"],
                        "atom_name": target_atom["atom_name"],
                        "atom_type": _decode_atom_type(encoded_matrix[target]),
                        "residue_name": target_atom["residue_name"],
                        "residue_number": target_atom["residue_number"],
                        "chain_id": target_atom["chain_id"],
                        "insertion_code": target_atom["insertion_code"],
                    },
                }
            )

    warning = mapping_warning
    if not chosen and warning is None:
        warning = "No positive gradient-weighted edge signal passed the interpretation filters."
    return _build_atom_interpretation(
        node_scores=node_scores,
        encoded_matrix=encoded_matrix,
        atom_metadata=atom_metadata,
        method=method,
        include_identifiers=include_identifiers,
        top_k=top_k,
        warning=warning,
        top_edges=top_edges,
    )


def _last_conv2d(model):
    last_layer = None
    for layer in model.modules():
        if isinstance(layer, nn.Conv2d):
            last_layer = layer
    return last_layer


def _prediction_payload(probabilities, experiment_scores, interpretation):
    if not probabilities:
        return None
    return {
        "mean_score": float(np.mean(probabilities)),
        "std_score": float(np.std(probabilities)),
        "experiment_breakdown": experiment_scores,
        "interpretation": interpretation,
    }


class CholNetBackend:
    def __init__(
        self,
        gat_path,
        gcn_path,
        gnn_path,
        k_ensembles=50,
        enable_vina_docking=None,
        docking_config=None,
    ):
        self.gat_path = gat_path
        self.gcn_path = gcn_path
        self.gnn_path = gnn_path
        self.k = k_ensembles
        self._inference_lock = threading.RLock()

        if enable_vina_docking is None:
            enabled_value = os.environ.get("CHOLNET_ENABLE_VINA_DOCKING", "1")
            enable_vina_docking = enabled_value.strip().lower() not in {
                "0", "false", "no", "off",
            }
        self.enable_vina_docking = bool(enable_vina_docking)
        self.docking_config = docking_config
        if self.enable_vina_docking and self.docking_config is None:
            self.docking_config = DockingConfig.from_environment(
                os.path.dirname(os.path.abspath(__file__))
            )

        self.gat_models = self._load_gat_models()
        self.gcn_models = self._load_gcn_models()
        self.gnn_models = self._load_gnn_models()

    def _load_gat_models(self):
        print("Loading GAT Models...")
        all_experiments = []
        for exp in range(1, 6):
            exp_models = []
            path_prefix = f"{self.gat_path}/GATModels-5A_exp{exp}v2/Models"
            for index in range(1, self.k + 1):
                model_path = f"{path_prefix}/model_bin_{index}.pth"
                if os.path.exists(model_path):
                    model = GAT(in_channels=37, out_channels=32).to(device)
                    model.load_state_dict(torch.load(model_path, map_location=device))
                    model.eval()
                    exp_models.append(model)
            all_experiments.append(exp_models)
        return all_experiments

    def _load_gcn_models(self):
        print("Loading GCN Models...")
        all_experiments = []
        for exp in range(1, 6):
            exp_models = []
            path_prefix = f"{self.gcn_path}/GCN-5A_Exp{exp}/Models"
            for index in range(1, self.k + 1):
                model_path = f"{path_prefix}/model_bin_{index}.pth"
                if os.path.exists(model_path):
                    model = GCN(input_dim=37).to(device)
                    model.load_state_dict(torch.load(model_path, map_location=device))
                    model.eval()
                    exp_models.append(model)
            all_experiments.append(exp_models)
        return all_experiments

    def _load_gnn_models(self):
        print("Loading GNN Models...")
        all_experiments = []
        for exp in range(1, 6):
            exp_models = []
            path_prefix = f"{self.gnn_path}/GNN-5A_Exp{exp}/Models"
            for index in range(1, self.k + 1):
                model_path = f"{path_prefix}/model_bin_{index}.pth"
                if os.path.exists(model_path):
                    model = CNN2D(input_channels=1).to(device)
                    model = nn.DataParallel(model)
                    model.load_state_dict(torch.load(model_path, map_location=device))
                    model.eval()
                    exp_models.append(model)
            all_experiments.append(exp_models)
        return all_experiments

    def predict(self, pdb_file_path, include_interpretation_ids=True):
        # Shared model objects receive gradients during interpretation. Serializing
        # requests prevents hooks/gradients from different Flask workers from mixing.
        with self._inference_lock:
            return self._predict_locked(pdb_file_path, include_interpretation_ids)

    def _predict_locked(self, pdb_file_path, include_interpretation_ids):
        if not os.path.exists(pdb_file_path):
            return {"status": "error", "message": "File not found"}

        session_id = str(uuid.uuid4())
        session_output_dir = os.path.join("TempProcessing", session_id)
        os.makedirs(session_output_dir, exist_ok=True)

        try:
            prediction_pdb_path = pdb_file_path
            docking_metadata = None
            docking_pose_by_clr = {}

            if self.enable_vina_docking:
                docking_result = dock_cholesterol(
                    pdb_file_path,
                    os.path.join(session_output_dir, "docking"),
                    config=self.docking_config,
                )
                prediction_pdb_path = docking_result.complex_pdb
                docking_metadata = docking_result.api_metadata()
                docking_pose_by_clr = {
                    (pose.residue_number, pose.chain_id): pose.as_dict()
                    for pose in docking_result.assigned_poses
                }

            preprocessing_dir = os.path.join(session_output_dir, "preprocessing")
            os.makedirs(preprocessing_dir, exist_ok=True)
            create_graphs([prediction_pdb_path], preprocessing_dir)

            graph_by_clr = {}
            matrix_by_clr = {}
            generated_pdb_paths = []
            for root, _, files in os.walk(preprocessing_dir):
                for filename in files:
                    full_path = os.path.join(root, filename)
                    if filename.lower().endswith(".pdb"):
                        generated_pdb_paths.append(full_path)
                    if filename.endswith("_graphs.npy"):
                        clr_key = parse_clr_from_path(full_path)
                        if clr_key:
                            graph_by_clr[clr_key] = full_path
                    elif filename.endswith("_combined_matrix.npy"):
                        clr_key = parse_clr_from_path(full_path)
                        if clr_key:
                            matrix_by_clr[clr_key] = full_path

            clr_keys = sorted(
                set(graph_by_clr) | set(matrix_by_clr),
                key=lambda key: (key[0], key[1]),
            )
            if not clr_keys:
                return {
                    "status": "error",
                    "message": "Preprocessing failed. No graph files were generated. "
                    + (
                        "No protein atoms were found within 5 Angstrom of any "
                        "docked cholesterol pose."
                        if self.enable_vina_docking
                        else "Ensure a correctly formatted CLR ligand is present."
                    ),
                }

            results = []
            for clr_key in clr_keys:
                graph_path = graph_by_clr.get(clr_key)
                matrix_path = matrix_by_clr.get(clr_key)
                encoded_matrix = None
                atom_metadata = None
                mapping_warning = "Graph data was unavailable for PDB atom mapping."

                if graph_path:
                    payload = _load_graph_payload(graph_path)
                    encoded_matrix = np.asarray(payload["encoded_matrix"], dtype=np.float32)
                    atom_metadata, mapping_warning = _map_graph_nodes_to_pdb(
                        original_pdb_path=prediction_pdb_path,
                        clr_key=clr_key,
                        encoded_matrix=encoded_matrix,
                        generated_pdb_paths=generated_pdb_paths,
                    )

                bucket = {
                    "filename": os.path.basename(pdb_file_path),
                    "clr_residue_number": clr_key[0],
                    "clr_chain_id": clr_key[1],
                    "site_source": (
                        "vina_docked"
                        if clr_key in docking_pose_by_clr
                        else (
                            "uploaded_experimental"
                            if self.enable_vina_docking
                            else "uploaded"
                        )
                    ),
                    "GAT": None,
                    "GCN": None,
                    "GNN": None,
                }
                if clr_key in docking_pose_by_clr:
                    bucket["docking_pose"] = docking_pose_by_clr[clr_key]

                if graph_path and encoded_matrix is not None:
                    bucket["GAT"] = self._evaluate_gat(
                        graph_path,
                        encoded_matrix,
                        atom_metadata,
                        mapping_warning,
                        include_interpretation_ids,
                    )
                    bucket["GCN"] = self._evaluate_gcn(
                        graph_path,
                        encoded_matrix,
                        atom_metadata,
                        mapping_warning,
                        include_interpretation_ids,
                    )
                gnn_encoded_matrix = encoded_matrix
                gnn_atom_metadata = atom_metadata
                gnn_mapping_warning = mapping_warning
                if matrix_path and gnn_encoded_matrix is None:
                    # Preserve score-only GNN inference if preprocessing produced
                    # the matrix but not its companion graph dictionary.
                    matrix_for_shape = np.load(matrix_path)
                    real_count = int(
                        np.count_nonzero(
                            ~np.isclose(matrix_for_shape, 0.0, atol=1e-8).all(axis=1)
                        )
                    )
                    gnn_encoded_matrix = np.zeros(
                        (real_count, len(ATOM_TYPE_NAMES)),
                        dtype=np.float32,
                    )
                    gnn_atom_metadata = None
                    gnn_mapping_warning = (
                        "The companion graph dictionary was unavailable; the GNN "
                        "score is valid but atom/residue interpretation is unavailable."
                    )

                if matrix_path and gnn_encoded_matrix is not None:
                    bucket["GNN"] = self._evaluate_gnn(
                        matrix_path,
                        gnn_encoded_matrix,
                        gnn_atom_metadata,
                        gnn_mapping_warning,
                        include_interpretation_ids,
                    )

                results.append(bucket)

            response = {
                "status": "success",
                "interpretation_schema_version": 1,
                "docking": docking_metadata,
                "results": results,
            }
            if include_interpretation_ids:
                # The single-file viewer needs the exact structure used for
                # prediction.  With Vina enabled, this is the untouched upload
                # (including its original HETATM ligands and residue names) plus
                # the newly docked CLR poses.
                with open(
                    prediction_pdb_path,
                    "r",
                    encoding="utf-8",
                    errors="replace",
                ) as handle:
                    response["structure_pdb"] = handle.read()
                response["structure_source"] = (
                    "uploaded_with_vina_clr_poses"
                    if self.enable_vina_docking
                    else "uploaded"
                )
            return response

        except DockingError as exc:
            return {"status": "error", "message": f"Vina docking failed: {exc}"}
        except Exception as exc:
            return {"status": "error", "message": str(exc)}
        finally:
            robust_rmtree(session_output_dir)

    def _evaluate_gat(
        self,
        file_path,
        encoded_matrix,
        atom_metadata,
        mapping_warning,
        include_identifiers,
    ):
        try:
            graph = organize_graph_gat(file_path, label=0).to(device)
        except Exception:
            return None

        batch = torch.zeros(graph.x.size(0), dtype=torch.long, device=device)
        probabilities = []
        experiment_scores = {}
        edge_runs = []

        for exp_index, models in enumerate(self.gat_models):
            exp_probabilities = []
            for model in models:
                model.zero_grad(set_to_none=True)
                try:
                    with torch.enable_grad():
                        output, (edge_index, alpha) = model(
                            graph.x,
                            graph.edge_index,
                            graph.edge_attr,
                            batch,
                        )
                        logit = output.reshape(-1)[0]
                        probability = float(torch.sigmoid(logit.detach()).item())
                        exp_probabilities.append(probability)
                        probabilities.append(probability)

                        if alpha.requires_grad:
                            alpha.retain_grad()
                            logit.backward()
                            if alpha.grad is not None:
                                importance = (alpha * alpha.grad).detach().cpu().numpy()
                                edge_runs.append(
                                    (edge_index.detach().cpu().numpy(), importance)
                                )
                except Exception as exc:
                    print(f"Warning: GAT interpretation skipped for one model: {exc}")
                finally:
                    model.zero_grad(set_to_none=True)

            if exp_probabilities:
                experiment_scores[f"exp_{exp_index + 1}"] = float(np.mean(exp_probabilities))

        interpretation = _edge_based_interpretation(
            edge_runs=edge_runs,
            encoded_matrix=encoded_matrix,
            atom_metadata=atom_metadata,
            method="Attention x gradient",
            include_identifiers=include_identifiers,
            mapping_warning=mapping_warning,
        )
        return _prediction_payload(probabilities, experiment_scores, interpretation)

    def _evaluate_gcn(
        self,
        file_path,
        encoded_matrix,
        atom_metadata,
        mapping_warning,
        include_identifiers,
    ):
        try:
            base_graph = organize_graph_gcn(file_path, label=0).to(device)
        except Exception:
            return None

        probabilities = []
        experiment_scores = {}
        edge_runs = []

        for exp_index, models in enumerate(self.gcn_models):
            exp_probabilities = []
            for model in models:
                model.zero_grad(set_to_none=True)
                try:
                    edge_attr = base_graph.edge_attr.detach().clone().requires_grad_(True)
                    graph = Data(
                        x=base_graph.x,
                        edge_index=base_graph.edge_index,
                        edge_attr=edge_attr,
                    )
                    graph.batch = torch.zeros(
                        graph.x.size(0), dtype=torch.long, device=device
                    )
                    with torch.enable_grad():
                        output = model(graph)
                        logit = output.reshape(-1)[0]
                        probability = float(torch.sigmoid(logit.detach()).item())
                        exp_probabilities.append(probability)
                        probabilities.append(probability)
                        logit.backward()
                        if edge_attr.grad is not None:
                            importance = (edge_attr * edge_attr.grad).detach().cpu().numpy()
                            edge_runs.append(
                                (
                                    graph.edge_index.detach().cpu().numpy(),
                                    importance,
                                )
                            )
                except Exception as exc:
                    print(f"Warning: GCN interpretation skipped for one model: {exc}")
                finally:
                    model.zero_grad(set_to_none=True)

            if exp_probabilities:
                experiment_scores[f"exp_{exp_index + 1}"] = float(np.mean(exp_probabilities))

        interpretation = _edge_based_interpretation(
            edge_runs=edge_runs,
            encoded_matrix=encoded_matrix,
            atom_metadata=atom_metadata,
            method="Edge weight x gradient",
            include_identifiers=include_identifiers,
            mapping_warning=mapping_warning,
        )
        return _prediction_payload(probabilities, experiment_scores, interpretation)

    def _evaluate_gnn(
        self,
        file_path,
        encoded_matrix,
        atom_metadata,
        mapping_warning,
        include_identifiers,
    ):
        try:
            matrix = np.load(file_path)
            if not isinstance(matrix, np.ndarray) or matrix.ndim != 2:
                return None
        except Exception:
            return None

        input_tensor = (
            torch.tensor(matrix, dtype=torch.float32)
            .unsqueeze(0)
            .unsqueeze(0)
            .to(device)
        )
        probabilities = []
        experiment_scores = {}
        cams = []

        for exp_index, models in enumerate(self.gnn_models):
            exp_probabilities = []
            for wrapped_model in models:
                model = wrapped_model.module if isinstance(wrapped_model, nn.DataParallel) else wrapped_model
                target_layer = _last_conv2d(model)
                activations = []
                gradients = []
                handles = []
                model.zero_grad(set_to_none=True)

                if target_layer is not None:
                    handles.append(
                        target_layer.register_forward_hook(
                            lambda _module, _inputs, output: activations.append(output)
                        )
                    )
                    handles.append(
                        target_layer.register_full_backward_hook(
                            lambda _module, _grad_input, grad_output: gradients.append(grad_output[0])
                        )
                    )

                try:
                    with torch.enable_grad():
                        output = model(input_tensor)
                        logit = output.reshape(-1)[0]
                        probability = float(torch.sigmoid(logit.detach()).item())
                        exp_probabilities.append(probability)
                        probabilities.append(probability)
                        logit.backward()

                        if activations and gradients:
                            activation = activations[-1]
                            gradient = gradients[-1]
                            weights = gradient.mean(dim=(2, 3), keepdim=True)
                            cam = F.relu((weights * activation).sum(dim=1, keepdim=True))
                            if cam.shape[-2:] != input_tensor.shape[-2:]:
                                cam = F.interpolate(
                                    cam,
                                    size=input_tensor.shape[-2:],
                                    mode="bilinear",
                                    align_corners=False,
                                )
                            cam_np = cam[0, 0].detach().cpu().numpy()
                            cams.append(_positive_normalize(cam_np))
                except Exception as exc:
                    print(f"Warning: GNN Grad-CAM skipped for one model: {exc}")
                finally:
                    for handle in handles:
                        handle.remove()
                    model.zero_grad(set_to_none=True)

            if exp_probabilities:
                experiment_scores[f"exp_{exp_index + 1}"] = float(np.mean(exp_probabilities))

        if cams:
            mean_cam = np.mean(np.stack(cams, axis=0), axis=0)
            node_scores = mean_cam[: encoded_matrix.shape[0]].sum(axis=1)
            warning = mapping_warning
        else:
            node_scores = np.zeros(encoded_matrix.shape[0], dtype=np.float64)
            warning = "Grad-CAM did not produce an importance map for this site."

        interpretation = _build_atom_interpretation(
            node_scores=node_scores,
            encoded_matrix=encoded_matrix,
            atom_metadata=atom_metadata,
            method="Grad-CAM",
            include_identifiers=include_identifiers,
            warning=warning,
        )
        return _prediction_payload(probabilities, experiment_scores, interpretation)