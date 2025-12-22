import glob
import os
import shutil
import uuid
import time
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.nn import GATConv, GCNConv, global_mean_pool
from collections import Counter

# Importing your custom preprocessing module
from FilterAtomsGraphs import create_graphs

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def robust_rmtree(path, retries=5, delay=0.2):
    if not os.path.exists(path):
        return

    for i in range(retries):
        try:
            shutil.rmtree(path)
            return
        except OSError:
            time.sleep(delay)

    print(f"Warning: Could not fully delete temp dir {path} due to file locks.")

class GAT(torch.nn.Module):
    def __init__(self, in_channels, out_channels, dropout_p=0.1):
        super().__init__()
        self.gat = GATConv(in_channels, out_channels, heads=1, concat=True, edge_dim=1)
        self.pool = global_mean_pool
        self.dropout = nn.Dropout(p=dropout_p)
        self.norm = nn.BatchNorm1d(out_channels)
        self.linear = torch.nn.Linear(out_channels, 1)

    def forward(self, x, edge_index, edge_attr, batch):
        out, attn_weights = self.gat(x, edge_index, edge_attr, return_attention_weights=True)
        out = self.dropout(out)
        out = self.pool(out, batch)
        out = self.norm(out)
        out = self.dropout(out)
        out = self.linear(out)
        return out, attn_weights


class GCN(nn.Module):
    def __init__(self, input_dim):
        super(GCN, self).__init__()
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
        x, edge_index, edge_weight, batch = data.x, data.edge_index, data.edge_attr, data.batch
        x = self.conv1(x, edge_index, edge_weight)
        x = self.bn1(x)
        x = F.relu(x)
        x = self.dropout_gcn(x)
        x = self.conv2(x, edge_index, edge_weight)
        x = self.bn2(x)
        x = F.relu(x)
        x = self.dropout_gcn(x)
        x = self.conv3(x, edge_index, edge_weight)
        x = self.bn3(x)
        x = F.relu(x)
        x = self.dropout_gcn(x)
        x = global_mean_pool(x, batch)
        x = self.fc1(x)
        x = F.relu(x)
        x = self.dropout(x)
        x = self.out(x)
        return x


class CNN2D(nn.Module):
    def __init__(self, input_channels):
        super(CNN2D, self).__init__()
        self.conv1 = nn.Sequential(
            nn.Conv2d(input_channels, 32, kernel_size=3, padding=1),
            nn.BatchNorm2d(32),
            nn.ReLU()
        )
        self.pool1 = nn.MaxPool2d(2, 2)
        self.conv2 = nn.Sequential(
            nn.Conv2d(32, 64, kernel_size=3, padding=1),
            nn.BatchNorm2d(64),
            nn.ReLU()
        )
        self.pool2 = nn.MaxPool2d(2, 2)
        self.conv3 = nn.Sequential(
            nn.Conv2d(64, 128, kernel_size=3, padding=1),
            nn.BatchNorm2d(128),
            nn.ReLU()
        )
        self.pool3 = nn.MaxPool2d(2, 2)
        self.flatten = nn.Flatten()
        self.fc1 = nn.Linear(128 * 4 * 18, 128)
        self.dropout = nn.Dropout(0.5)
        self.fc2 = nn.Linear(128, 1)

    def forward(self, x):
        x = self.conv1(x)
        x = self.pool1(x)
        x = self.conv2(x)
        x = self.pool2(x)
        x = self.conv3(x)
        x = self.pool3(x)
        x = self.flatten(x)
        x = torch.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.fc2(x)
        return x

def organize_graph_gat(file_path, label=0):
    data = np.load(file_path, allow_pickle=True).item()
    inverse_distance = data['inverse_distance']
    encoded_matrix = data['encoded_matrix']
    x = torch.tensor(encoded_matrix, dtype=torch.float32)
    adj = torch.tensor(inverse_distance, dtype=torch.float32)
    adj = adj / (adj.sum(dim=1, keepdim=True) + 1e-8)
    edge_index = (adj > 0).nonzero(as_tuple=False).t()
    edge_weight = adj[adj > 0]
    y = torch.tensor([label], dtype=torch.float32)
    return Data(x=x, edge_index=edge_index, edge_attr=edge_weight, y=y)


def organize_graph_gcn(file_path, label=0):
    data = np.load(file_path, allow_pickle=True).item()
    inverse_distance = data['inverse_distance']
    encoded_matrix = data['encoded_matrix']
    x = torch.tensor(encoded_matrix, dtype=torch.float32)
    adj = torch.tensor(inverse_distance, dtype=torch.float32)
    edge_index = (adj > 0).nonzero(as_tuple=False).t()
    edge_weight = adj[adj > 0]
    y = torch.tensor([label], dtype=torch.float32)
    return Data(x=x, edge_index=edge_index, edge_attr=edge_weight, y=y)

class CholNetBackend:
    def __init__(self, gat_path, gcn_path, gnn_path, k_ensembles=50):
        self.gat_path = gat_path
        self.gcn_path = gcn_path
        self.gnn_path = gnn_path
        self.k = k_ensembles

        self.gat_models = self._load_gat_models()
        self.gcn_models = self._load_gcn_models()
        self.gnn_models = self._load_gnn_models()

    def _load_gat_models(self):
        print("Loading GAT Models...")
        all_experiments = []
        for exp in range(1, 6):
            exp_models = []
            path_prefix = f"{self.gat_path}/GATModels-5A_exp{exp}v2/Models"
            for i in range(1, self.k + 1):
                model_path = f"{path_prefix}/model_bin_{i}.pth"
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
            for i in range(1, self.k + 1):
                model_path = f"{path_prefix}/model_bin_{i}.pth"
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
            for i in range(1, self.k + 1):
                model_path = f"{path_prefix}/model_bin_{i}.pth"
                if os.path.exists(model_path):
                    model = CNN2D(input_channels=1).to(device)
                    model = nn.DataParallel(model)
                    model.load_state_dict(torch.load(model_path, map_location=device))
                    model.eval()
                    exp_models.append(model)
            all_experiments.append(exp_models)
        return all_experiments

    def predict(self, pdb_file_path):
        if not os.path.exists(pdb_file_path):
            return {"status": "error", "message": "File not found"}

        session_id = str(uuid.uuid4())
        session_output_dir = os.path.join("TempProcessing", session_id)
        os.makedirs(session_output_dir, exist_ok=True)

        try:
            create_graphs([pdb_file_path], session_output_dir)

            npy_files = []
            for root, dirs, files in os.walk(session_output_dir):
                for file in files:
                    if file.endswith(".npy"):
                        npy_files.append(os.path.join(root, file))

            if not npy_files:
                return {"status": "error",
                        "message": "Preprocessing failed. No graph files generated. Ensure CLR ligand is present."}

            merged_results = {
                "filename": os.path.basename(pdb_file_path),
                "GAT": None,
                "GCN": None,
                "GNN": None
            }

            for file_path in npy_files:
                gat = self._evaluate_gat(file_path)
                if gat:
                    merged_results["GAT"] = gat

                gcn = self._evaluate_gcn(file_path)
                if gcn:
                    merged_results["GCN"] = gcn

                gnn = self._evaluate_gnn(file_path)
                if gnn:
                    merged_results["GNN"] = gnn

            return {"status": "success", "results": [merged_results]}

        except Exception as e:
            return {"status": "error", "message": str(e)}

        finally:
            # Cleanup using the robust method to handle Windows locks
            robust_rmtree(session_output_dir)

    def _evaluate_gat(self, file_path):
        try:
            graph = organize_graph_gat(file_path, label=0).to(device)
        except:
            return None

        total_probs = []
        experiment_scores = {}

        for exp_idx, models in enumerate(self.gat_models):
            exp_probs = []
            for model in models:
                with torch.no_grad():
                    out, _ = model(
                        graph.x,
                        graph.edge_index,
                        graph.edge_attr,
                        batch=torch.zeros(graph.x.size(0), dtype=torch.long).to(device)
                    )
                    exp_probs.append(torch.sigmoid(out).item())

            avg_exp = np.mean(exp_probs)
            experiment_scores[f"exp_{exp_idx + 1}"] = float(avg_exp)
            total_probs.extend(exp_probs)

        return {
            "mean_score": float(np.mean(total_probs)),
            "std_score": float(np.std(total_probs)),
            "experiment_breakdown": experiment_scores
        }

    def _evaluate_gcn(self, file_path):
        try:
            graph = organize_graph_gcn(file_path, label=0).to(device)
            graph.batch = torch.zeros(graph.num_nodes, dtype=torch.long).to(device)
        except:
            return None

        total_probs = []
        experiment_scores = {}

        for exp_idx, models in enumerate(self.gcn_models):
            exp_probs = []
            for model in models:
                with torch.no_grad():
                    out = model(graph)
                    exp_probs.append(torch.sigmoid(out).item())

            avg_exp = np.mean(exp_probs)
            experiment_scores[f"exp_{exp_idx + 1}"] = float(avg_exp)
            total_probs.extend(exp_probs)

        return {
            "mean_score": float(np.mean(total_probs)),
            "std_score": float(np.std(total_probs)),
            "experiment_breakdown": experiment_scores
        }

    def _evaluate_gnn(self, file_path):
        try:
            grid = np.load(file_path)
            if isinstance(grid, dict) or grid.ndim != 2:
                return None
        except:
            return None

        grid_tensor = torch.tensor(grid, dtype=torch.float32).unsqueeze(0).unsqueeze(0)
        grid_tensor = grid_tensor.to(device)

        total_probs = []
        experiment_scores = {}

        for exp_idx, models in enumerate(self.gnn_models):
            exp_probs = []
            for model in models:
                with torch.no_grad():
                    output = model(grid_tensor).squeeze(1)
                    exp_probs.append(torch.sigmoid(output).item())

            avg_exp = np.mean(exp_probs)
            experiment_scores[f"exp_{exp_idx + 1}"] = float(avg_exp)
            total_probs.extend(exp_probs)

        return {
            "mean_score": float(np.mean(total_probs)),
            "std_score": float(np.std(total_probs)),
            "experiment_breakdown": experiment_scores
        }