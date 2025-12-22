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
import EnergyCode
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

    if not file.filename.endswith('.pdb'):
        return jsonify({"status": "error", "message": "Invalid file type. Please upload a .pdb file."}), 400

    session_id = str(uuid.uuid4())
    session_dir = os.path.join(UPLOAD_FOLDER, session_id)
    os.makedirs(session_dir, exist_ok=True)

    filepath = os.path.join(session_dir, file.filename)

    try:
        file.save(filepath)

        response = cholnet_backend.predict(filepath)
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


if __name__ == '__main__':
    #app.run(debug=True)
    # port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=5000, debug=True, use_reloader=False)