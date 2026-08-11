import sys
print(sys.executable)
from flask import Flask, request, jsonify, render_template, send_from_directory
from flask_cors import CORS
import numpy as np
import copy
import os
from itertools import islice
import matplotlib
#matplotlib.use('Agg')  # Use a non-interactive backend
import matplotlib.pyplot as plt
from collections import Counter
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import Graph2D
import SubtractedCorrelationMatrix
import SaltBridgePlot
import json
import logging
import uuid
import time
import shutil
import PDBCompareMethods
import Allosteric
from cholbind import CholNetBackend  # same import you used before

app=Flask(__name__)
application=app
CORS(app, resources={r"/api/*": {"origins": ["http://flowallostery.westernu.edu", "http://204.9.174.78", "http://localhost:5173"]}})
logging.basicConfig(filename='/home/alexhernandez/Allosteric-Network/server/logfile.log', level=logging.ERROR)

#387,388,389,389,390,391,392
#328,329,334,338,378,348

UPLOAD_FOLDER = 'uploads'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

def robust_rmtree(path, retries=5, delay=0.2):
    if not os.path.exists(path):
        return
    for _ in range(retries):
        try:
            shutil.rmtree(path)
            return
        except OSError:
            time.sleep(delay)
    print(f"Warning: Could not delete upload session {path}")

print("Initializing CholNet Backend...")
cholnet_backend = CholNetBackend(
    gat_path="CholBindModels/GAT_Model",
    gcn_path="CholBindModels/GCN_Model",
    gnn_path="CholBindModels/GNN_Model",
    k_ensembles=5
)
print("CholNet Backend Initialized and Ready.")

@app.route('/api/py3dmol', methods=['POST'])
def py3dmol():
    #py3dmol_content = EnergyCode.visualizeBetweenness()

    with open('view_data.json', 'r') as json_file:
        py3dmol_content = json.load(json_file)

    return jsonify(py3dmol_content)

@app.route('/api/upload', methods=['POST'])
def upload_file():
    graph_data = Graph2D.process_graph_data()
    return jsonify(graph_data)

@app.route('/api/subtract', methods=['POST'])
def subtracted_matrix():
    plots_and_protein = SubtractedCorrelationMatrix.get_plots_and_protein_structure()
    return jsonify(plots_and_protein)

@app.route('/api/rerender', methods=['POST'])
def rerender():
    structure = SubtractedCorrelationMatrix.recalculate_from_new_cutoff_value()
    return jsonify(structure)

@app.route('/api/rerender-salt', methods=['POST'])
def rerender_salt():
    structure = SaltBridgePlot.recalculate_from_new_cutoff_value()
    return jsonify(structure)

@app.route('/api/extract_chains', methods=['POST'])
def extract_chains_from_pdb():
    return PDBCompareMethods.extract_chains()

@app.route('/api/allosteric', methods=['POST'])
def allosteric_plots():
    return Allosteric.process_graph_data()

@app.route('/api/calculate', methods=['POST'])
def compute_flow_betweenness():
    global adj_matrix
    global source_array
    global sink_array

    resist_1 = int(request.form['resist1'])
    resist_2 = int(request.form['resist2'])
    
    tempAdjDense = adj_matrix.todense()
    
    # Convert adjacency matrix to Laplacian matrix
    tempLapDense = -copy.deepcopy(tempAdjDense)
    for ii, irow in enumerate(tempLapDense):
        tempLapDense[ii, ii] = -np.sum(irow)
    
    # Compute pseudoinverse of Laplacian matrix
    tempLinv = np.linalg.pinv(tempLapDense)
    
    total_betweenness_score = 0

    for s in source_array:
        for t in sink_array:
            # Compute flow betweenness for the given edge
            v_source_sink_resist1 = tempLinv[resist_1, s] - tempLinv[resist_1, t]
            v_source_sink_resist2 = tempLinv[resist_2, s] - tempLinv[resist_2, t]
            b_resist1_resist2 = tempAdjDense[resist_1, resist_2] * (v_source_sink_resist1 - v_source_sink_resist2)
            total_betweenness_score += b_resist1_resist2

    #divide by number of combinations
    num_of_combinations = len(source_array) * len(sink_array)
    total_betweenness_score /= num_of_combinations

    betweenness_score = total_betweenness_score.item() # Convert to a standard Python type
    if betweenness_score < 0:
        betweenness_score *= -1
    
    return jsonify({'betweenness_score': betweenness_score})

@app.route('/api/cholnet', methods=['POST'])
def cholnet_predict():
    if 'pdb_file' not in request.files:
        return jsonify({"status": "error", "message": "No file part: pdb_file"}), 400

    file = request.files['pdb_file']
    if file.filename == '':
        return jsonify({"status": "error", "message": "No selected file"}), 400

    if not file.filename.lower().endswith('.pdb'):
        return jsonify({"status": "error", "message": "Invalid file type. Please upload a .pdb file."}), 400

    session_id = str(uuid.uuid4())
    session_dir = os.path.join(UPLOAD_FOLDER, session_id)
    os.makedirs(session_dir, exist_ok=True)

    safe_filename = os.path.basename(file.filename)
    filepath = os.path.join(session_dir, safe_filename)

    try:
        file.save(filepath)

        response = cholnet_backend.predict(
            filepath,
            include_interpretation_ids=True,
        )
        # Expecting your backend returns:
        # {"status":"success","results":[...]} or {"status":"error","message":"..."}
        if response.get("status") == "success":
            return jsonify(response), 200
        else:
            return jsonify({
                "status": "error",
                "message": response.get("message", "Processing error.")
            }), 500

    except Exception as e:
        app.logger.exception("CholNet predict failed")
        return jsonify({"status": "error", "message": f"Server Error: {str(e)}"}), 500

    finally:
        robust_rmtree(session_dir)

def _summarize_batch_interpretations(items, top_k=10):
    """Aggregate general residue/atom types without exposing structure IDs."""
    model_names = ("GNN", "GAT", "GCN")
    accumulators = {
        model_name: {
            "method": None,
            "sites_interpreted": 0,
            "atom_types": {},
            "residue_types": {},
        }
        for model_name in model_names
    }

    def add_records(target, records, label_key):
        for record in records or []:
            label = record.get(label_key)
            if not label:
                continue
            count = int(record.get("count", 0) or 0)
            if count <= 0:
                continue
            importance = float(record.get("mean_importance", 0.0) or 0.0)
            bucket = target.setdefault(
                label,
                {"count": 0, "weighted_importance": 0.0},
            )
            bucket["count"] += count
            bucket["weighted_importance"] += importance * count

    for item in items:
        for site in item.get("results", []):
            for model_name in model_names:
                model_result = site.get(model_name) or {}
                interpretation = model_result.get("interpretation") or {}
                if not interpretation:
                    continue
                accumulator = accumulators[model_name]
                accumulator["sites_interpreted"] += 1
                accumulator["method"] = interpretation.get("method")
                add_records(
                    accumulator["atom_types"],
                    interpretation.get("top_atom_types"),
                    "atom_type",
                )
                add_records(
                    accumulator["residue_types"],
                    interpretation.get("top_residue_types"),
                    "residue_name",
                )

    def finalize(values, label_key):
        total = sum(value["count"] for value in values.values())
        rows = []
        for label, value in values.items():
            count = value["count"]
            rows.append(
                {
                    label_key: label,
                    "count": count,
                    "percent": (count / total * 100.0) if total else 0.0,
                    "mean_importance": value["weighted_importance"] / count,
                }
            )
        rows.sort(
            key=lambda row: (
                -row["count"],
                -row["mean_importance"],
                row[label_key],
            )
        )
        return rows[:top_k]

    models = {}
    for model_name, accumulator in accumulators.items():
        models[model_name] = {
            "method": accumulator["method"],
            "sites_interpreted": accumulator["sites_interpreted"],
            "top_atom_types": finalize(accumulator["atom_types"], "atom_type"),
            "top_residue_types": finalize(
                accumulator["residue_types"],
                "residue_name",
            ),
        }

    return {
        "mode": "batch_general",
        "files_processed": sum(1 for item in items if item.get("results")),
        "models": models,
    }


def _remove_site_interpretation_details(items):
    """Batch mode returns only the aggregate type summary, never PDB/node IDs."""
    for item in items:
        for site in item.get("results", []):
            for model_name in ("GNN", "GAT", "GCN"):
                model_result = site.get(model_name)
                if isinstance(model_result, dict):
                    model_result.pop("interpretation", None)

@app.route('/api/cholnet/batch', methods=['POST'])
def cholnet_predict_batch():
    # Must match frontend: formData.append("pdb_files", f)
    files = request.files.getlist('pdb_files')
    if not files or len(files) == 0:
        return jsonify({"status": "error", "message": "No files part: pdb_files"}), 400

    session_id = str(uuid.uuid4())
    session_dir = os.path.join(UPLOAD_FOLDER, session_id)
    os.makedirs(session_dir, exist_ok=True)

    items = []
    try:
        for f in files:
            # Some browsers include directory in filename when using webkitdirectory
            # e.g. "subdir/foo.pdb" -> keep only basename
            original_name = f.filename or ""
            basename = os.path.basename(original_name)

            if basename == "":
                items.append({"filename": original_name, "results": [], "error": "Empty filename"})
                continue

            if not basename.lower().endswith(".pdb"):
                items.append({"filename": basename, "results": [], "error": "Invalid file type (not .pdb)"})
                continue

            filepath = os.path.join(session_dir, basename)

            try:
                f.save(filepath)

                response = cholnet_backend.predict(
                    filepath,
                    include_interpretation_ids=False,
                )
                if response.get("status") == "success":
                    items.append({
                        "filename": basename,
                        "results": response.get("results", [])
                    })
                else:
                    items.append({
                        "filename": basename,
                        "results": [],
                        "error": response.get("message", "Processing error.")
                    })

            except Exception as e:
                app.logger.exception(f"CholNet batch predict failed for {basename}")
                items.append({
                    "filename": basename,
                    "results": [],
                    "error": f"Server Error: {str(e)}"
                })

        interpretation_summary = _summarize_batch_interpretations(items)
        _remove_site_interpretation_details(items)

        return jsonify({
            "status": "success",
            "items": items,
            "interpretation_summary": interpretation_summary,
        }), 200

    except Exception as e:
        app.logger.exception("CholNet batch predict failed (outer)")
        return jsonify({"status": "error", "message": f"Server Error: {str(e)}"}), 500

    finally:
        robust_rmtree(session_dir)


if __name__ == '__main__':
    #app.run(debug=True)
    # port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=5000, debug=True, use_reloader=False)