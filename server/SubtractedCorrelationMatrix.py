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

# index for residue dictionary
RESIDUE_NAME = 0
X = 1
Y = 2
Z = 3

def calculate_cb_distance_matrix(topology_file):
    """
    Calculates a distance matrix based on the Cβ atoms from an MD trajectory.
    """
    # Load the MD trajectory
    u = mda.Universe(topology_file)

    # Select Cβ atoms (excluding Glycine)
    cb_selection = u.select_atoms("name CB and not resname GLY")
    
    # # If you want to include Cα for Glycine
    # ca_selection = u.select_atoms("name CA and resname GLY")
    
    # Combine Cβ and Cα for Glycine
    all_cb_atoms = cb_selection 

    # Initialize an empty distance matrix
    num_atoms = len(all_cb_atoms)
    distance_matrix = np.zeros((num_atoms, num_atoms))

    # Loop through each frame in the trajectory to calculate distances
    for ts in u.trajectory:
        # Get current positions of Cβ atoms
        cb_positions = all_cb_atoms.positions
        
        # Calculate pairwise distances
        for i in range(num_atoms):
            for j in range(num_atoms):
                distance_matrix[i, j] += np.linalg.norm(cb_positions[i] - cb_positions[j])
    
    # Average over all frames
    distance_matrix /= len(u.trajectory)
    count = np.sum(~np.isnan(distance_matrix) & (distance_matrix != 0))
    print(count, " is count")
    
    return distance_matrix

# Assuming we have two matrices of different shapes and we want to perform element-wise subtraction

def subtract_overlapping_regions(matrix1, matrix2):
    """
    Performs element-wise subtraction of overlapping regions between two matrices.
    If the matrices are of different shapes, only the overlapping region is subtracted.
    """
    # Determine the shape of the overlapping region
    min_rows = min(matrix1.shape[0], matrix2.shape[0]) 
    print ("min_rows:", min_rows)
    min_cols = min(matrix1.shape[1], matrix2.shape[1])
    print ("min_cols:",min_cols)
    
    # Perform element-wise subtraction on the overlapping region
    result = matrix1[:min_rows, :min_cols] - matrix2[:min_rows, :min_cols]
    
    return result

def create_edgelist_from_mda_universe(residue_pairs, residue_dictionary):
    global RESIDUE_NAME, X, Y, Z

    edge_list = []

    for resID1, resID2 in residue_pairs:
        residue1 = residue_dictionary.get(resID1)
        residue2 = residue_dictionary.get(resID2)
        
        # Use MDAnalysis to calculate the center of mass for each residue
        crd1 = {residue1[X], residue1[Y], residue1[Z]}
        crd2 = {residue2[X], residue2[Y], residue2[Z]}
        
        resname1 = residue1[RESIDUE_NAME]
        resid1 = resID1
        
        resname2 = residue2[RESIDUE_NAME]
        resid2 = resID2
        
        # Create an edge label based on the residue names and IDs
        edgeLabel = f'{resname1}.{resid1}-{resname2}.{resid2} ({resID1-1}-{resID2-1})'

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
    return edge_list[:10000]
    
def create_residue_pairs_list(csv_file):
    # Load the CSV file
    df = pd.read_csv(csv_file)
    
    # Extract the relevant columns: ResidueID1 and ResidueID2
    residue_pairs = df[['ResidueID1', 'ResidueID2']].values.tolist()
    
    return residue_pairs

def create_res_dictionary(csv_file):
    # Create an empty dictionary to store Residue ID as key and (Residue Name, X, Y, Z) as values
    residue_dict = {}

    # Read data from the CSV file
    with open(csv_file, mode='r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header row

        for row in reader:
            # Extract relevant information from each row
            residue_id = int(row[3])  # Residue ID
            residue_name = row[2]      # Residue Name
            x = float(row[5])          # X coordinate
            y = float(row[6])          # Y coordinate
            z = float(row[7])          # Z coordinate

            # Store the values in the dictionary
            residue_dict[residue_id] = (residue_name, x, y, z)
    return residue_dict


def get_plots():
    pdb_file1 = request.files['pdb_file1']
    pdb_file2 = request.files['pdb_file2']

    pdb_file1_path = 'pdb_file1.pdb'
    pdb_file2_path = 'pdb_file2.pdb'
    pdb_file1.save(pdb_file1_path)
    pdb_file2.save(pdb_file2_path)

    cb_distance_matrixA = calculate_cb_distance_matrix(pdb_file1_path)

    # residue_pairs = create_residue_pairs_list("test.csv")
    # residue_dictionary = create_res_dictionary("WT_CB_distance_pairs.csv")
    # edge_list = create_edgelist_from_mda_universe(residue_pairs, residue_dictionary)

    ## checking purpose
    print (np.shape(cb_distance_matrixA))

    cb_distance_matrixB=calculate_cb_distance_matrix(pdb_file2_path)

    # Example matrices (replace these with actual matrices from the notebook code)
    matrix1 = cb_distance_matrixA
    matrix2 = cb_distance_matrixB

    # Apply the element-wise subtraction function to the matrices
    result = np.abs(subtract_overlapping_regions(matrix1, matrix2))

    # Example matrices (replace these with actual matrices from the notebook code)
    matrix1 = cb_distance_matrixA
    matrix2 = cb_distance_matrixB

    # Apply the element-wise subtraction function to the matrices
    result = np.abs(subtract_overlapping_regions(matrix1, matrix2))

    # Assuming matrix1, matrix2, and result are your data matrices

    # Create a figure with three subplots (1 row, 3 columns)
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))  # Adjust figsize for a better layout

    # Plot matrix1
    cax1 = axs[0].imshow(matrix1, cmap="viridis", interpolation="nearest")
    fig.colorbar(cax1, ax=axs[0], label="Distance (Å)")
    axs[0].set_title("WT")
    axs[0].set_xlabel("Residue Index")
    axs[0].set_ylabel("Residue Index")
    axs[0].invert_yaxis()  # Invert y-axis as before

    # Plot matrix2
    cax2 = axs[1].imshow(matrix2, cmap="viridis", interpolation="nearest")
    fig.colorbar(cax2, ax=axs[1], label="Distance (Å)")
    axs[1].set_title("Mut")
    axs[1].set_xlabel("Residue Index")
    axs[1].invert_yaxis()  # Invert y-axis as before

    # Plot the thresholded matrix (result)
    cax3 = axs[2].imshow(result, cmap="viridis", interpolation="nearest")
    fig.colorbar(cax3, ax=axs[2], label="Distance (Å)")
    axs[2].set_title("Subtracted Wt and Mut")
    axs[2].set_xlabel("Residue Index")
    axs[2].invert_yaxis()  # Invert y-axis as before

    # Adjust layout for better spacing between plots
    plt.tight_layout()

    buffer = io.BytesIO()
    plt.savefig(buffer, format="png")
    buffer.seek(0)
    calculated_matrix_image = base64.b64encode(buffer.getvalue()).decode('utf8')
    buffer.close()
    plt.close()

    # Assuming 'sub' is your subtracted matrix
    threshold = 3

    # Set all values below or equal to 5 to NaN (or 0 if you prefer)
    sub_thresholded = np.where(result > threshold, result, np.nan)  # Use np.nan for better visual distinction in the heatmap

    print(sub_thresholded)

    num_residues = result.shape[0]
    print ("number residue:",num_residues)
    # Plot the thresholded matrix
    plt.imshow(sub_thresholded, cmap="viridis", interpolation="nearest")
    plt.colorbar(label="Distance (Å)")
    plt.title("Subtracted Distance Matrix (Threshold > 5Å)")
    plt.xlabel("Residue Index")
    plt.ylabel("Residue Index")
    # Invert y-axis as before
    plt.gca().invert_yaxis()

    buffer = io.BytesIO()
    plt.savefig(buffer, format="png")
    buffer.seek(0)
    subtracted_distance_matrix_image = base64.b64encode(buffer.getvalue()).decode('utf8')
    buffer.close()
    plt.close()

    with open(pdb_file2_path, 'r') as file:
        pdb_content = file.read()

    # view_data = {
    #     'file_content': file_content,
    #     'edges': edge_list
    # }

    plots = {
        'calculated_matrix_image' : calculated_matrix_image,
        'subtracted_distance_matrix_image' : subtracted_distance_matrix_image,
        'pdb_content' : pdb_content
        # 'edges' : edge_list
    }

    return plots