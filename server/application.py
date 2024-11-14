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

app=Flask(__name__)
application=app
CORS(app, resources={r"/api/*": {"origins": ["http://flowallostery.westernu.edu", "http://204.9.174.78", "http://localhost:5173"]}})
logging.basicConfig(filename='/home/alexhernandez/Allosteric-Network/server/logfile.log', level=logging.ERROR)

#387,388,389,389,390,391,392
#328,329,334,338,378,348

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

if __name__ == '__main__':
    #app.run(debug=True)
    # port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=5000, debug=True)
