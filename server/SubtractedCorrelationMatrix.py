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

#Add original atom index to be able to do new matrix
def pdb_to_dataframe(pdb_file):
    """
    Load a PDB file using MDAnalysis and convert key atom information to a pandas DataFrame.
    """
    u = mda.Universe(pdb_file)
    
    # Extract atom-related data: atom name, residue name, residue ID, and chain ID
    atom_data = {
        'Atom Name': u.atoms.names,
        'Residue Name': u.atoms.resnames,
        'Residue ID': u.atoms.resids,
        'Chain ID': u.atoms.segids,
        'X': u.atoms.positions[:, 0],
        'Y': u.atoms.positions[:, 1],
        'Z': u.atoms.positions[:, 2],
    }
    
    # Create a pandas DataFrame from the atom data
    df = pd.DataFrame(atom_data)
    
    return df

def get_residue_ids(pdb_file):
    """
    Extracts the residue IDs from a given PDB file using MDAnalysis.
    """
    u = mda.Universe(pdb_file)
    residues = u.select_atoms("name CA or (name CB and not resname GLY)").residues
    residue_ids = residues.resids  # Extract the residue IDs
    return residue_ids

def filter_by_residue_ids(df, common_residue_ids):
    """
    Filters the DataFrame to include only rows with Residue IDs that are in the common_residue_ids list.
    """
    # Filter the DataFrame for rows where the 'Residue ID' is in the list of common residue IDs
    filtered_df = df[df['Residue ID'].isin(common_residue_ids)]
    
    return filtered_df

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

    # Calculate pairwise distances using broadcasting
    diff = coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]
    distances = np.sqrt(np.sum(diff ** 2, axis=-1))

    # Create a DataFrame for the distance matrix with proper indexing
    distance_df = pd.DataFrame(distances, index=df['NewIndex'], columns=df['NewIndex'])

    return distance_df

# Function to update the universe atom coordinates
def update_coordinates_in_universe(u, df, sub, hashmap_cb1):
    # Scale down shifts to prevent large movements
    shift_scale = 0.1

    # Go through each index in the distance matrix
    for new_index in sub.index:
        # Get the corresponding Residue ID and Chain ID using the hashmap
        residue_id, chain_id = hashmap_cb1[new_index]

        # Select all atoms for this residue in the MDAnalysis universe
        residue_atoms = u.select_atoms(f"resid {residue_id} and segid {chain_id}")
        
        # Ensure we have selected atoms
        if len(residue_atoms) > 0:
            # Calculate the shift for this residue based on the subtracted matrix
            shift_vector = np.sum(sub.loc[new_index].values) * shift_scale
            
            # Apply the shift equally to all atoms in this residue
            # This part assumes the shift vector affects the X, Y, Z directions equally
            residue_atoms.positions += np.array([shift_vector, shift_vector, shift_vector])

def recalculate_from_new_cutoff_value():
    pdb_file1 = request.files['pdb_file1']
    pdb_file2 = request.files['pdb_file2']

    energy_value = float(request.form['energy'])

    pdb_file1_path = 'pdb_file1.pdb'
    pdb_file2_path = 'pdb_file2.pdb'
    pdb_file1.save(pdb_file1_path)
    pdb_file2.save(pdb_file2_path)

    u = mda.Universe(pdb_file1_path)

    residue_pairs = create_residue_pairs_list("merged_distance_pairs.csv", energy_value)
    edge_list = create_edgelist_from_mda_universe_and_residue_pairs(u, residue_pairs)

    with open(pdb_file1_path, 'r') as file:
        pdb_content = file.read()

    # view_data = {
    #     'file_content': file_content,
    #     'edges': edge_list
    # }

    structure = {
        'pdb_content' : pdb_content,
        'edges' : edge_list
    }

    return structure

def create_edgelist_from_mda_universe_and_residue_pairs(pubStrucUniverse, residue_pairs):
    global RESIDUE_NAME, X, Y, Z

    edge_list = []

    for resID1, chainID1, resID2 , chainID2 in residue_pairs[:10]:
        residue1 = pubStrucUniverse.select_atoms(f"resid {resID1} and segid {chainID1}")
        residue2 = pubStrucUniverse.select_atoms(f"resid {resID2} and segid {chainID2}")

        empty_residue = False

        # Check if the atom groups are empty by their length
        if len(residue1) == 0 or len(residue2) == 0:
            empty_residue = True
            print("test if empty")

        if not empty_residue:
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

def residue_pairs_for_sub(hash, sub):
    # filtered_sub = sub.where((sub >= 3.0) & (sub <= 8.0), 0)
    # residue_pairs = []
    # for index, row in filtered_sub.iterrows():
    #     for col_index, value in row.items():
    #         print(f"Row {index}, Column {col_index}: {value}")
    #         if not math.isnan(value):
    #             resID1, chainID1 = hash.get(index)
    #             resID2, chainID2 = hash.get(col_index)
    #             residue_pairs.append(resID1, chainID1, resID2, chainID2)
    # return residue_pairs
    # Filter the sub DataFrame and stack it to get pairs of indices with non-zero values
    filtered_sub = sub[(sub >= 3.0) & (sub <= 8.0)].stack()
    
    # List comprehension to extract residue pairs using the hashmap
    residue_pairs = [
        (hash.get(index1) + hash.get(index2))  # Combine tuples to create a flat structure
        for (index1, index2), value in filtered_sub.items()
        if hash.get(index1) and hash.get(index2)  # Ensure both indices exist in hashmap
    ]
    print(len(residue_pairs), " is length of residue pairs")
    return residue_pairs

def create_residue_pairs_list(csv_file, distance = 6.0):
    # Load the CSV file
    df = pd.read_csv(csv_file)
    
    # Determine cutoff
    filtered_df = df.loc[df['Distance'] >= distance]

    residue_pairs = filtered_df[['ResidueID1', 'ChainID1', 'ResidueID2', 'ChainID2']].values.tolist()
    
    return residue_pairs

def get_plots(pdb_file1_path, pdb_file2_path):
# Load your first PDB file
    pdb_file1 = pdb_file1_path  
    sys1 = pdb_to_dataframe(pdb_file1)
    # sys1['System']='WT'
    sys1 = sys1.reset_index()
    print(sys1.head())  

    u = mda.Universe(pdb_file1)

    # Load your second PDB file 
    pdb_file2 = pdb_file2_path 
    sys2 = pdb_to_dataframe(pdb_file2)
    # sys2['System']='Mut'
    sys2=sys2.reset_index()
    print(sys2.tail())  

    # Extract residue IDs from both PDB files
    residue_ids_1 = get_residue_ids(pdb_file1)
    residue_ids_2 = get_residue_ids(pdb_file2)

    # Find the common residue IDs between both PDB files
    common_residue_ids = np.intersect1d(residue_ids_1, residue_ids_2)

    # Print or return the common residues
    print("Common Residue IDs between the two PDB files:", common_residue_ids)
    # output_file = "common_residues.dat"
    # np.savetxt(output_file, common_residue_ids, fmt='%d')

    # Filter df1 and df2 for common residue IDs
    # extract CB/GLY CA data for each residue
    filtered_df1 = filter_by_residue_ids(sys1, common_residue_ids)
    filtered_cb1 = filtered_df1.query('`Atom Name` == "CB" | (`Atom Name` == "CA" & `Residue Name` == "GLY")')


    filtered_df2 = filter_by_residue_ids(sys2, common_residue_ids)
    filtered_cb2 = filtered_df2.query('`Atom Name` == "CB" | (`Atom Name` == "CA" & `Residue Name` == "GLY")')
    # Reset the index and keep it as a column
    # Display the filtered DataFrames
    print("Filtered df1 with common residue IDs:")
    print(filtered_df1)

    print("\nFiltered df2 with common residue IDs:")
    print(filtered_df2)

    print(filtered_cb1['Residue ID'].head(), " is head")

    # re-index, RESIDUE ID and CHAIN ID ARE STILL THE SAME VALUES AFTER THIS
    filtered_cb1['NewIndex']=range(0,len(filtered_cb1))
    filtered_cb2['NewIndex']=range(0,len(filtered_cb2))

    print(len(filtered_cb1), " is length of f_cb1")

    print(filtered_cb1['Residue ID'].head(), " is head. Test")

    # Create hashmaps for filtered_cb1 and filtered_cb2 to return Residue ID and Chain ID from NewIndex
    hashmap_cb1 = {row['NewIndex']: (row['Residue ID'], row['Chain ID']) for _, row in filtered_cb1.iterrows()}
    hashmap_cb2 = {row['NewIndex']: (row['Residue ID'], row['Chain ID']) for _, row in filtered_cb2.iterrows()}

    matrixA=compute_pairwise_distances(filtered_cb1)
    print(matrixA.head(), " is matrix a head")
    matrixB=compute_pairwise_distances(filtered_cb2)
    sub=np.abs(matrixA-matrixB)

    # Apply the updates based on the subtracted distances
    update_coordinates_in_universe(u, filtered_cb1, sub, hashmap_cb1) 
    u.atoms.write("updated_structure.pdb")

    new_universe = mda.Universe("updated_structure.pdb")  

    #get residue pairs for edgelist
    residue_pairs = residue_pairs_for_sub(hashmap_cb1, sub)
    #new_edgelist = create_edgelist_from_mda_universe_and_residue_pairs(u, residue_pairs)

    print(sub.head(), " is sub head")

    # Create a figure with three subplots (1 row, 3 columns)
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))  # Adjust figsize for a better layout

    # Plot matrix1
    cax1 = axs[0].imshow(matrixA, cmap="viridis", interpolation="nearest")
    fig.colorbar(cax1, ax=axs[0], label="Distance (Å)")
    axs[0].set_title("WT")
    axs[0].set_xlabel("Residue Index")
    axs[0].set_ylabel("Residue Index")
    axs[0].invert_yaxis()  # Invert y-axis as before

    # Plot matrix2
    cax2 = axs[1].imshow(matrixB, cmap="viridis", interpolation="nearest")
    fig.colorbar(cax2, ax=axs[1], label="Distance (Å)")
    axs[1].set_title("Mut")
    axs[1].set_xlabel("Residue Index")
    axs[1].invert_yaxis()  # Invert y-axis as before

    # Plot the thresholded matrix (result)
    cax3 = axs[2].imshow(sub, cmap="viridis", interpolation="nearest")
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

    # save WT and Mut CB distance paris to a dataframe table
    filtered_cb1.to_csv('WT_CB_distance_pairs.csv',index=False)
    filtered_cb2.to_csv('Mut_CB_distance_pairs.csv',index=False)

    # filter out the distance between 5 to 10A from subtracted matrix
    threshold_min = 3
    threshold_max = 8

    # Set all values below or equal to 5 to NaN (or 0 if you prefer)
    sub_thresholded = np.where((sub > threshold_min) & (sub < threshold_max), sub, np.nan)  # Use np.nan for better visual distinction in the heatmap

    print(sub_thresholded, " is thresh")

    num_residues = sub.shape[0]
    print ("number residue:",num_residues)
    # Plot the thresholded matrix
    plt.imshow(sub_thresholded, cmap="viridis", interpolation="nearest")
    plt.colorbar(label="Distance (Å)")
    plt.title("Subtracted Distance Matrix (3<x<8A)")
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

    ## show the top 5 largest distance from subtracted matrix ## 
    residue_ids = np.arange(1, sub.shape[0])  # Assuming residue IDs are 1 to N

    # Convert to a DataFrame if needed
    df = pd.DataFrame(sub, index=residue_ids, columns=residue_ids)

    flattened_matrix = df.unstack().reset_index()  
    flattened_matrix.columns = ['Residue1', 'Residue2', 'Distance']

    # Keep only unique pairs (Residue1 < Residue2) to avoid duplicates 
    unique_pairs = flattened_matrix[flattened_matrix['Residue1'] < flattened_matrix['Residue2']]

    # filter distance within limitation 
    filter_distance=unique_pairs.query('Distance >3 & Distance <8')
    print(filter_distance.head(), " is head")
    print(filter_distance.tail(), " is tail")
    #print (filter_distance)
    # Sort by distance to get the top 5 largest values
    top_5 = filter_distance.nlargest(5, 'Distance')

    # Extract distances as a 1D array
    distances = filter_distance['Distance'].values

    # Calculate the number of bins using the square root choice
    num_bins = int(np.ceil(np.sqrt(len(distances))))

    # Plot the distribution of filtered distances
    plt.hist(distances, bins=num_bins, edgecolor='black')
    plt.xlabel('Distance')
    plt.ylabel('Frequency')
    plt.title(f'Distribution of Filtered Distances (3 < Distance < 8)')

    buffer = io.BytesIO()
    plt.savefig(buffer, format="png")
    buffer.seek(0)
    distribution_graph = base64.b64encode(buffer.getvalue()).decode('utf8')
    buffer.close()
    plt.close()

    # Show the result
    print(top_5)

    #save to a dataframe table
    top_5.to_csv('Top5_distance_pairs.csv',index=False)

    return calculated_matrix_image, subtracted_distance_matrix_image, distribution_graph#, new_edgelist


def get_plots_and_protein_structure():
    pdb_file1 = request.files['pdb_file1']
    pdb_file2 = request.files['pdb_file2']

    pdb_file1_path = 'pdb_file1.pdb'
    pdb_file2_path = 'pdb_file2.pdb'
    pdb_file1.save(pdb_file1_path)
    pdb_file2.save(pdb_file2_path)

    calculated_matrix_image, subtracted_distance_matrix_image, distribution_graph = get_plots(pdb_file1_path, pdb_file2_path)

    u = mda.Universe(pdb_file1_path)
    print(len(u.atoms), " is length of atoms")

    residue_pairs = create_residue_pairs_list("merged_distance_pairs.csv")
    print(residue_pairs[:10], " are first 10 resi`due pairs from csv")
    edge_list = create_edgelist_from_mda_universe_and_residue_pairs(u, residue_pairs)

    with open(pdb_file1_path, 'r') as file:
        pdb_content = file.read()

    # view_data = {
    #     'file_content': file_content,
    #     'edges': edge_list
    # }

    plots = {
        'calculated_matrix_image' : calculated_matrix_image,
        'subtracted_distance_matrix_image' : subtracted_distance_matrix_image,
        'distribution_graph' : distribution_graph,
        'pdb_content' : pdb_content,
        'edges' : edge_list
    }

    return plots