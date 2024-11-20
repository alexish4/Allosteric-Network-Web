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
import math
import json
import SaltBridgePlot
import PDBCompareMethods

def get_residue_ids(pdb_file):
    """
    Extracts the residue IDs from a given PDB file using MDAnalysis.
    """
    u = mda.Universe(pdb_file)
    residues = u.select_atoms("name CA or (name CB and not resname GLY)").residues
    residue_ids = residues.resids  # Extract the residue IDs
    return residue_ids

def compute_pairwise_distances(df):
    """
    Compute the pairwise Euclidean distances between residues based on their 3D coordinates.

    Parameters:
    df (pd.DataFrame): DataFrame containing 'X', 'Y', 'Z' coordinates and 'NewIndex' as index.

    Returns:
    pd.DataFrame: A DataFrame containing the pairwise distance matrix.
    """
    # Extract the coordinates (X, Y, Z)
    coordinates = df[['X', 'Y', 'Z']].values
    residue_ids = df['Residue ID'].values
    chain_ids = df['Chain ID'].values
    new_index = df['NewIndex'].values

    # Calculate pairwise distances using broadcasting
    diff = coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]
    distances = np.sqrt(np.sum(diff ** 2, axis=-1))

    # Create a DataFrame for the distance matrix with proper indexing
    #distance_df = pd.DataFrame(distances, index=df['NewIndex'], columns=df['NewIndex'])

    # Initialize an empty list to store the results
    results = []

    # Iterate over the pairs of residues and chains
    for i in range(len(df)):
        for j in range(i + 1, len(df)):
            # Collect the necessary information
            residue_id1 = residue_ids[i]
            residue_id2 = residue_ids[j]
            chain_id1 = chain_ids[i]
            chain_id2 = chain_ids[j]
            distance = distances[i, j]
            new_index1 = new_index[i]
            new_index2 = new_index[j]

            # Append the information as a tuple
            results.append([residue_id1, residue_id2, chain_id1, chain_id2, distance, new_index1, new_index2])

    # Create a DataFrame from the results
    distance_df = pd.DataFrame(results, columns=['ResidueID1', 'ResidueID2', 'ChainID1', 'ChainID2', 'Distance', 'Index1', 'Index2'])

    return distance_df
    

def recalculate_from_new_cutoff_value():
    pdb_file1 = request.files['pdb_file1']
    pdb_file2 = request.files['pdb_file2']
    edge_file = request.files.get('edge_file')

    # Check if a file was submitted
    submitted_edge_file = False
    if edge_file and edge_file.filename != '':
        submitted_edge_file = True
        edge_file.save("Subtract_Files/edges_table.csv")
    
    file_to_render = "Subtract_Files/saved_sub1.csv"
    if submitted_edge_file:
        filter_by_edge_file("Subtract_Files/saved_sub1.csv", "Subtract_Files/edges_table.csv")
        file_to_render = "Subtract_Files/filtered_edges.csv"

    lower_bound = float(request.form['lower_bound'])

    selected_chains = json.loads(request.form.get('selected_chains'))
    chain_ranges = json.loads(request.form.get('chain_ranges'))

    validated_ranges = {}
    for chain, ranges_str in chain_ranges.items():
        if ranges_str:
            validated_ranges[chain] = PDBCompareMethods.parse_ranges(ranges_str)
    print(validated_ranges)

    filtered_chains = [chain for chain, is_selected in selected_chains.items() if is_selected]

    pdb_file1_path = 'Subtract_Files/pdb_file1.pdb'
    pdb_file2_path = 'Subtract_Files/pdb_file2.pdb'
    pdb_file1.save(pdb_file1_path)
    pdb_file2.save(pdb_file2_path)

    u = mda.Universe(pdb_file1_path)

    residue_pairs, subtracted_table = PDBCompareMethods.create_residue_pairs_list(file_to_render, filtered_chains, validated_ranges, lower_bound)
    print(len(residue_pairs), " is length of residue pairs")
    edge_list = PDBCompareMethods.rerender_edgelist_from_mda_universe_and_residue_pairs(u, residue_pairs)

    with open(pdb_file1_path, 'r') as file:
        pdb_content = file.read()

    print(subtracted_table.columns, " is subtracted table columns")

    structure = {
        'pdb_content' : pdb_content,
        'edges' : edge_list,
        'table' : subtracted_table.to_json(orient='records') 
    }

    return structure

def filter_by_edge_file(current_csv, csv_with_edges_to_filter_by):
    current_df = pd.read_csv(current_csv)
    edges_df = pd.read_csv(csv_with_edges_to_filter_by)

    # Making sure format matches
    first_chain_id = current_df['ChainID1'].iloc[0]
    if first_chain_id in {'PROA', 'PROB', 'PROC', 'PROD'}:
        edges_df['ChainID1'] = edges_df['ChainID1'].replace({'A': 'PROA', 'B': 'PROB', 'C': 'PROC', 'D': 'PROD'})
        edges_df['ChainID2'] = edges_df['ChainID2'].replace({'A': 'PROA', 'B': 'PROB', 'C': 'PROC', 'D': 'PROD'})

    edges_df = edges_df.drop("Distance", axis=1, errors="ignore")

    # Merging the two DataFrames on the specified columns
    filtered_df = current_df.merge(
        edges_df, 
        on=["ResidueID1", "ChainID1", "ResidueID2", "ChainID2"],
        how="inner"
    )

    # Selecting only the columns from current_df
    filtered_df = filtered_df[current_df.columns]
    filtered_df.to_csv("Subtract_Files/filtered_edges.csv", index=False)
    
    return filtered_df

def get_plots(pdb_file1_path, pdb_file2_path):
# Load your first PDB file
    pdb_file1 = pdb_file1_path  
    sys1 = PDBCompareMethods.pdb_to_dataframe(pdb_file1)
    # sys1['System']='WT'
    sys1 = sys1.reset_index()

    # Load your second PDB file 
    pdb_file2 = pdb_file2_path 
    sys2 = PDBCompareMethods.pdb_to_dataframe(pdb_file2)
    sys2=sys2.reset_index()

    # Extract residue IDs from both PDB files
    residue_ids_1 = get_residue_ids(pdb_file1)
    residue_ids_2 = get_residue_ids(pdb_file2)

    # Find the common residue IDs between both PDB files
    common_residue_ids = np.intersect1d(residue_ids_1, residue_ids_2)

    # Filter df1 and df2 for common residue IDs
    # extract CB/GLY CA data for each residue
    filtered_df1 = PDBCompareMethods.filter_by_residue_ids(sys1, common_residue_ids)
    filtered_cb1 = filtered_df1.query('`Atom Name` == "CB" | (`Atom Name` == "CA" & `Residue Name` == "GLY")')

    filtered_df2 = PDBCompareMethods.filter_by_residue_ids(sys2, common_residue_ids)
    filtered_cb2 = filtered_df2.query('`Atom Name` == "CB" | (`Atom Name` == "CA" & `Residue Name` == "GLY")')

    # re-index, RESIDUE ID and CHAIN ID ARE STILL THE SAME VALUES AFTER THIS
    filtered_cb1['NewIndex']=range(0,len(filtered_cb1))
    filtered_cb2['NewIndex']=range(0,len(filtered_cb2))

    matrixA=compute_pairwise_distances(filtered_cb1)
    matrixB=compute_pairwise_distances(filtered_cb2)

    sub = {'Index1': matrixA['Index1'],
    'Index2':matrixA['Index2'],
    'Distance_wt':matrixA['Distance'],
    'Distance_mut':matrixB['Distance'],
    'Delta_Distance':np.abs(matrixA['Distance']-matrixB['Distance']),
    'ResidueID1':matrixA['ResidueID1'],
    'ChainID1':matrixA['ChainID1'],
    'ResidueID2':matrixA['ResidueID2'],
    'ChainID2':matrixA['ChainID2']}
    sub=pd.DataFrame(sub)

    filtered_sub = sub[(sub['Distance_wt'] < 15) & (sub['Distance_mut'] < 15)]

    filtered_sub.to_csv("Subtract_Files/saved_sub1.csv", index=False)

    # Create a figure with three subplots (1 row, 3 columns)
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))  # Adjust figsize for a better layout

    # Plot matrix1
    sc1 = axs[0].imshow(matrixA.pivot('Index1', 'Index2', 'Distance'), cmap='viridis', aspect='auto')
    cbar1 = fig.colorbar(sc1, ax=axs[0])
    cbar1.set_label('Distance (Å)', fontsize=14)
    axs[0].set_title("PDB 1", fontsize=16)
    axs[0].set_xlabel("Residue Index", fontsize=14)
    axs[0].set_ylabel("Residue Index", fontsize=14)
    axs[0].invert_yaxis()

    # Plot matrix2
    sc2 = axs[1].imshow(matrixB.pivot('Index1', 'Index2', 'Distance'), cmap='viridis', aspect='auto')
    cbar2 = fig.colorbar(sc2, ax=axs[1])
    cbar2.set_label('Distance (Å)', fontsize=14)
    axs[1].set_title("PDB 2", fontsize=16)
    axs[1].set_xlabel("Residue Index", fontsize=14)
    axs[1].invert_yaxis()

    # Plot the thresholded matrix (result)
    sc3 = axs[2].imshow(sub.pivot('Index1', 'Index2', 'Delta_Distance'), cmap='viridis', aspect='auto')
    cbar3 = fig.colorbar(sc3, ax=axs[2])
    cbar3.set_label('∆ Distance (Å)', fontsize=14)
    axs[2].set_title("∆ Distance", fontsize=16)
    axs[2].set_xlabel("Residue Index", fontsize=14)
    axs[2].invert_yaxis()

    # Adjust layout for better spacing between plots
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.4)

    buffer = io.BytesIO()
    plt.savefig(buffer, format="png")
    buffer.seek(0)
    calculated_matrix_image = base64.b64encode(buffer.getvalue()).decode('utf8')
    buffer.close()
    plt.close()

    # Extract distances as a 1D array
    #distances = matrixA.values
    distances = filtered_sub['Delta_Distance'].values

    # Calculate the number of bins using the square root choice
    num_bins = int(np.ceil(np.sqrt(len(distances))))

    # Plot the distribution of filtered distances
    plt.hist(distances, bins=num_bins, edgecolor='black')
    plt.xlabel('∆ Distance (Å)', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.title(f'Distribution of ∆ Distances', fontsize=16)

    buffer = io.BytesIO()
    plt.savefig(buffer, format="png")
    buffer.seek(0)
    distribution_graph = base64.b64encode(buffer.getvalue()).decode('utf8')
    buffer.close()
    plt.close()

    return calculated_matrix_image, distribution_graph


def get_plots_and_protein_structure():
    pdb_file1 = request.files['pdb_file1']
    pdb_file2 = request.files['pdb_file2']

    pdb_file1_path = 'Subtract_Files/pdb_file1.pdb'
    pdb_file2_path = 'Subtract_Files/pdb_file2.pdb'
    pdb_file1.save(pdb_file1_path)
    pdb_file2.save(pdb_file2_path)

    calculated_matrix_image, distribution_graph = get_plots(pdb_file1_path, pdb_file2_path)
    salt_plot_image, salt_distribution_image = SaltBridgePlot.generate_salt_plot(pdb_file1_path, pdb_file2_path)

    with open(pdb_file1_path, 'r') as file:
        pdb_content = file.read()

    edge_list = []

    plots = {
        'calculated_matrix_image' : calculated_matrix_image,
        'salt_plot_image' : salt_plot_image,
        'subtract_distribution_graph' : distribution_graph,
        'salt_distribution_graph' : salt_distribution_image,
        'pdb_content' : pdb_content,
        'edges' : edge_list
    }

    return plots