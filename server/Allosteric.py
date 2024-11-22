from flask import request, jsonify
import uuid

def get_plots_and_protein_structure():
    pdb_file = request.files['pdb_file']

    # Generate unique filenames
    unique_id = uuid.uuid4().hex  # Generate a unique identifier
    pdb_file_path = f'Subtract_Files/{unique_id}_pdb_file1.pdb'
    pdb_file.save(pdb_file_path)

    with open(pdb_file_path, 'r') as file:
        pdb_content = file.read()

    edge_list = []

    plots = {
        'pdb_content' : pdb_content,
        'edges' : edge_list,
        'unique_id' : unique_id
    }

    return jsonify(plots)