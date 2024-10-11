import MDAnalysis as mda
from MDAnalysis.analysis.distances import self_distance_array
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform, pdist  # Import squareform from scipy
from bokeh.plotting import figure, show
from bokeh.models import HoverTool, ColorBar, LinearColorMapper
from bokeh.io import output_notebook
from bokeh.transform import linear_cmap
from bokeh.models.tickers import BasicTicker
from bokeh.models.formatters import PrintfTickFormatter
import numpy as np
import matplotlib.pyplot as plt
import io
from flask import request
from Bio.PDB import PDBParser, Select
import base64
import pandas as pd
import csv

def create_edgelist_from_mda_universe_and_residue_pairs(pubStrucUniverse, residue_pairs):
    global RESIDUE_NAME, X, Y, Z

    edge_list = []

    for resID1, chainID1, resID2 , chainID2 in residue_pairs:
        residue1 = pubStrucUniverse.select_atoms(f"resid {resID1} and segid {chainID1}")
        residue2 = pubStrucUniverse.select_atoms(f"resid {resID2} and segid {chainID2}")

        # Use MDAnalysis to calculate the center of mass for each residue
        crd1 = residue1.center_of_mass()
        crd2 = residue2.center_of_mass()

        resname1 = residue1.residues[0].resname
        resid1 = int(residue1.residues[0].resid)
        
        resname2 = residue2.residues[0].resname
        resid2 = int(residue2.residues[0].resid)
        
        # Create an edge label based on the residue names and IDs
        edgeLabel = f'{resname1}.{resid1}-{resname2}.{resid2} ({resID1}.{chainID1}-{resID2}.{chainID2})'

        #converting to python types instead of numpy types so we can jsonify
        edge_data = {
            'label': edgeLabel,
            'coords': {
                'start': [float(c) for c in crd1],  # Convert NumPy array to Python list of floats
                'end': [float(c) for c in crd2]  # Convert NumPy array to Python list of floats
            }
        }

        edge_list.append(edge_data)
    print(len(edge_list), " is length of edge list")
    return edge_list
    
def create_residue_pairs_list(csv_file):
    # Load the CSV file
    df = pd.read_csv(csv_file)
    
    # Determine cutoff
    filtered_df = df.loc[df['Distance'] >= 6.0]

    residue_pairs = filtered_df[['ResidueID1', 'ChainID1', 'ResidueID2', 'ChainID2']].values.tolist()
    
    return residue_pairs

def get_plots():
    pdb_file1 = request.files['pdb_file1']
    pdb_file2 = request.files['pdb_file2']

    pdb_file1_path = 'pdb_file1.pdb'
    pdb_file2_path = 'pdb_file2.pdb'
    pdb_file1.save(pdb_file1_path)
    pdb_file2.save(pdb_file2_path)

    u = mda.Universe(pdb_file1_path)

    residue_pairs = create_residue_pairs_list("merged_distance_pairs.csv")
    edge_list = create_edgelist_from_mda_universe_and_residue_pairs(u, residue_pairs)

    with open(pdb_file1_path, 'r') as file:
        pdb_content = file.read()

    # view_data = {
    #     'file_content': file_content,
    #     'edges': edge_list
    # }

    plots = {
        # 'calculated_matrix_image' : calculated_matrix_image,
        # 'subtracted_distance_matrix_image' : subtracted_distance_matrix_image,
        'pdb_content' : pdb_content,
        'edges' : edge_list
    }

    return plots