import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform  # Import squareform from scipy

import numpy as np

import matplotlib.pyplot as plt
import pandas as pd
import io
import base64
import SubtractedCorrelationMatrix
from flask import request
import json
import PDBCompareMethods

def get_residue_ids(pdb_file):
    """
    Extracts the residue IDs from a given PDB file using MDAnalysis.
    """
    u = mda.Universe(pdb_file)
    #residues = u.select_atoms("name CA or (name CB and not resname GLY)").residues
    residues = u.select_atoms(
        "(name CZ and resname ARG) or "
        "(name NZ and resname LYS) or "
        "(name CG and resname ASP) or "
        "(name CD and resname GLU)"
    ).residues

    residue_ids = residues.resids  # Extract the residue IDs
    return residue_ids

def recalculate_from_new_cutoff_value():
    pdb_file1 = request.files['pdb_file1']
    pdb_file2 = request.files['pdb_file2']
    edge_file = request.files.get('edge_file')

    # Check if a file was submitted
    submitted_edge_file = False
    if edge_file and edge_file.filename != '':
        submitted_edge_file = True
        edge_file.save("Subtract_Files/edges_table.csv")
    
    file_to_render = "Subtract_Files/salt_bridge.csv"
    if submitted_edge_file:
        SubtractedCorrelationMatrix.filter_by_edge_file("Subtract_Files/saved_sub.csv", "Subtract_Files/edges_table.csv")
        file_to_render = "Subtract_Files/filtered_edges.csv"

    lower_bound = float(request.form['lower_bound'])

    selected_chains = json.loads(request.form.get('selected_chains'))
    chain_ranges = json.loads(request.form.get('chain_ranges'))

    validated_ranges = {}
    for chain, ranges_str in chain_ranges.items():
        if ranges_str:
            validated_ranges[chain] = PDBCompareMethods.parse_ranges(ranges_str)
    print(validated_ranges)

    filtered_chains = []
    if selected_chains.get('A'):
        filtered_chains.append("A")
        filtered_chains.append("PROA")
    if selected_chains.get('B'):
        filtered_chains.append("B")
        filtered_chains.append("PROB")
    if selected_chains.get('C'):
        filtered_chains.append("C")
        filtered_chains.append("PROC")
    if selected_chains.get('D'):
        filtered_chains.append("D")
        filtered_chains.append("PROD")

    pdb_file1_path = 'Subtract_Files/pdb_file1.pdb'
    pdb_file2_path = 'Subtract_Files/pdb_file2.pdb'
    pdb_file1.save(pdb_file1_path)
    pdb_file2.save(pdb_file2_path)

    u = mda.Universe(pdb_file1_path)

    residue_pairs, salt_table = PDBCompareMethods.create_residue_pairs_list(file_to_render, filtered_chains, validated_ranges, lower_bound)
    print(len(residue_pairs), " is length of residue pairs")
    edge_list = PDBCompareMethods.rerender_edgelist_from_mda_universe_and_residue_pairs(u, residue_pairs)

    with open(pdb_file1_path, 'r') as file:
        pdb_content = file.read()

    print(salt_table.columns, " is salt table columns")

    structure = {
        'pdb_content' : pdb_content,
        'edges' : edge_list,
        'table' : salt_table.to_json(orient='records') 
    }

    return structure

def generate_salt_plot(pdb_file1, pdb_file2):
    sys1 = PDBCompareMethods.pdb_to_dataframe(pdb_file1)
    sys1 = sys1.reset_index()

    sys2 = PDBCompareMethods.pdb_to_dataframe(pdb_file2)
    sys2 = sys2.reset_index()

    u_wt=mda.Universe(pdb_file1, pdb_file1)
    u1=mda.Universe(pdb_file2, pdb_file2)

    # Extract residue IDs from both PDB files
    residue_ids_1 = get_residue_ids(pdb_file1)
    residue_ids_2 = get_residue_ids(pdb_file2)

    # Find the common residue IDs between both PDB files
    common_residue_ids = np.intersect1d(residue_ids_1, residue_ids_2)

    filtered_df1 = PDBCompareMethods.filter_by_residue_ids(sys1, common_residue_ids)

    filtered_cb1 = filtered_df1.query(
        '(`Atom Name` == "CZ" & `Residue Name` == "ARG") | '
        '(`Atom Name` == "NZ" & `Residue Name` == "LYS") | '
        '(`Atom Name` == "CG" & `Residue Name` == "ASP") | '
        '(`Atom Name` == "CD" & `Residue Name` == "GLU")'
    )


    filtered_df2 = PDBCompareMethods.filter_by_residue_ids(sys2, common_residue_ids)
    filtered_cb2 = filtered_df2.query(
        '(`Atom Name` == "CZ" & `Residue Name` == "ARG") | '
        '(`Atom Name` == "NZ" & `Residue Name` == "LYS") | '
        '(`Atom Name` == "CG" & `Residue Name` == "ASP") | '
        '(`Atom Name` == "CD" & `Residue Name` == "GLU")'
    )

    filtered_cb1['NewIndex']=range(1,len(filtered_cb1)+1)
    filtered_cb2['NewIndex']=range(1,len(filtered_cb2)+1)

    selection_string = " or ".join([
        f"(resid {res} and ((name CZ and resname ARG) or (name NZ and resname LYS) or "
        f"(name CG and resname ASP) or (name CD and resname GLU)))"
        for res in common_residue_ids
    ])
    # Select the CB atoms (or CA for GLY residues) from the specified residues
    atoms = u_wt.select_atoms(selection_string)

    # Container to store pairwise distances for each frame
    all_distances = []

    # Calculate pairwise distances across all frames
    for ts in u_wt.trajectory:
        # Get the positions of selected atoms for the current frame
        coordinates = atoms.positions
        #print (coordinates.shape)
        # Calculate pairwise distances (condensed distance array)
        distance_matrix = distances.self_distance_array(coordinates)
        
        # Convert the condensed distance array to a square matrix
        distance_square = squareform(distance_matrix)
        
        # Create a DataFrame with proper indexing for easy interpretation
        distance_df = pd.DataFrame(distance_square, index=atoms.indices, columns=atoms.indices)
        
        # Append the DataFrame for the current frame to the list
        all_distances.append(distance_df)
    
    selection_string = " or ".join([
        f"(resid {res} and ((name CZ and resname ARG) or (name NZ and resname LYS) or "
        f"(name CG and resname ASP) or (name CD and resname GLU)))"
        for res in common_residue_ids
    ])

    # Select the CB atoms (or CA for GLY residues) from the specified residues
    atoms_mut = u1.select_atoms(selection_string)

    # Container to store pairwise distances for each frame
    all_distances_mut = []

    # Calculate pairwise distances across all frames
    for ts_mut in u1.trajectory:
        # Get the positions of selected atoms for the current frame
        coordinates = atoms_mut.positions
        #print (coordinates.shape)
        # Calculate pairwise distances (condensed distance array)
        distance_matrix = distances.self_distance_array(coordinates)
        
        # Convert the condensed distance array to a square matrix
        distance_square = squareform(distance_matrix)
        
        # Create a DataFrame with proper indexing for easy interpretation
        distance_df = pd.DataFrame(distance_square, index=atoms_mut.indices, columns=atoms_mut.indices)
        
        # Append the DataFrame for the current frame to the list
        all_distances_mut.append(distance_df)

    # Convert to a numpy array for easy averaging
    all_distances = np.array(all_distances)
    # Calculate the average distance across frames
    avg_dist_arr = np.mean(all_distances, axis=0)
    std_dist_arr=np.std(all_distances,axis=0)

    # Convert to a numpy array for easy averaging
    all_distances_mut = np.array(all_distances_mut)
    # Calculate the average distance across frames
    avg_dist_arr_mut = np.mean(all_distances_mut, axis=0)
    std_dist_arr_mut=np.std(all_distances_mut,axis=0)

    ## show the top 5 largest distance from subtracted matrix ## 
    residue_ids = np.arange(1, avg_dist_arr.shape[0]+1)  # Assuming residue IDs are 1 to N

    # Convert to a DataFrame if needed
    df = pd.DataFrame(avg_dist_arr, index=residue_ids, columns=residue_ids)
    df_std=pd.DataFrame(std_dist_arr, index=residue_ids,columns=residue_ids)
    flattened_matrix = df.unstack().reset_index()  
    flattened_matrix.columns = ['Index1', 'Index2', 'Distance']

    flattened_matrix_std = df_std.unstack().reset_index()  
    flattened_matrix_std.columns = ['Index1', 'Index2', 'std']
    # Keep only unique pairs (Residue1 < Residue2) to avoid duplicates 
    unique_pairs = flattened_matrix[flattened_matrix['Index1'] < flattened_matrix['Index2']]
    unique_pairs_std = flattened_matrix_std[flattened_matrix_std['Index1'] < flattened_matrix_std['Index2']]
    unique_pairs['std']=unique_pairs_std['std']

    residue_ids_mut = np.arange(1, avg_dist_arr_mut.shape[0]+1)  # Assuming residue IDs are 1 to N

    # Convert to a DataFrame if needed
    df_mut = pd.DataFrame(avg_dist_arr_mut, index=residue_ids_mut, columns=residue_ids_mut)
    df_std_mut=pd.DataFrame(std_dist_arr_mut, index=residue_ids_mut,columns=residue_ids_mut)
    flattened_matrix_mut = df_mut.unstack().reset_index()  
    flattened_matrix_mut.columns = ['Index1', 'Index2', 'Distance']

    flattened_matrix_std_mut = df_std_mut.unstack().reset_index()  
    flattened_matrix_std_mut.columns = ['Index1', 'Index2', 'std']
    # Keep only unique pairs (Residue1 < Residue2) to avoid duplicates
    unique_pairs_mut = flattened_matrix_mut[flattened_matrix_mut['Index1'] < flattened_matrix_mut['Index2']]
    unique_pairs_std_mut = flattened_matrix_std_mut[flattened_matrix_std_mut['Index1'] < flattened_matrix_std_mut['Index2']]
    unique_pairs_mut['std']=unique_pairs_std_mut['std']

    merged_table=pd.merge(unique_pairs, filtered_cb1, left_on='Index1', right_on='NewIndex', how='inner')

    merged_mut=pd.merge(unique_pairs_mut, filtered_cb2, left_on='Index1', right_on='NewIndex', how='inner')

    merged_table.drop(['index','Atom Name','X','Y','Z','NewIndex'],axis=1,inplace=True)
    merged_mut.drop(['index','Atom Name','X','Y','Z','NewIndex'],axis=1,inplace=True)

    merged_table.rename(columns={'Residue Name':'ResidueName1',
                        'Residue ID':'ResidueID1',
                        'Chain ID': 'ChainID1'},inplace=True)
    merged_mut.rename(columns={'Residue Name':'ResidueName1',
                        'Residue ID':'ResidueID1',
                        'Chain ID': 'ChainID1'},inplace=True)
    merged_b=pd.merge(merged_table, filtered_cb1, left_on='Index2', right_on='NewIndex', how='inner')
    merged_bmut=pd.merge(merged_mut, filtered_cb2, left_on='Index2', right_on='NewIndex', how='inner')

    merged_b.drop(['index','Atom Name','NewIndex','X','Y','Z'],axis=1,inplace=True)
    merged_b.rename(columns={'Residue Name':'ResidueName2',
                        'Residue ID':'ResidueID2',
                        'Chain ID': 'ChainID2'
                            },inplace=True)
    merged_bmut.drop(['index','Atom Name','NewIndex','X','Y','Z'],axis=1,inplace=True)
    merged_bmut.rename(columns={'Residue Name':'ResidueName2',
                        'Residue ID':'ResidueID2',
                        'Chain ID': 'ChainID2'
                            },inplace=True)
    
    sub = {'Index1': merged_b['Index1'],
        'Index2':merged_b['Index2'],
        'Distance_wt':merged_b['Distance'],
        'Distance_mut':merged_bmut['Distance'],
        'Delta_Distance':np.abs(merged_b['Distance']-merged_bmut['Distance']),
        'ResidueName1':merged_b['ResidueName1'],
        'ResidueID1':merged_b['ResidueID1'],
        'ChainID1':merged_b['ChainID1'],
        'ResidueName2':merged_b['ResidueName2'],
        'ResidueID2':merged_b['ResidueID2'],
        'ChainID2':merged_b['ChainID2'],
        'STD_WT': merged_b['std'],
        'STD_MUT': merged_bmut['std']}
    sub=pd.DataFrame(sub)

    fig, axs = plt.subplots(1, 3, figsize=(18, 6))  # Adjust figsize for a better layout

    sc1 = axs[0].scatter(merged_b['Index1'], merged_b['Index2'], c=merged_b['Distance'], cmap='viridis')
    cbar1 = fig.colorbar(sc1, ax=axs[0])
    cbar1.set_label('Z Value')
    axs[0].set_xlabel('Index ID')
    axs[0].set_ylabel('Index ID')
    axs[0].set_title('WT-MUT Matrix Plot')

    # Second scatter plot
    sc2 = axs[1].scatter(merged_bmut['Index1'], merged_bmut['Index2'], c=merged_bmut['Distance'], cmap='viridis')
    cbar2 = fig.colorbar(sc2, ax=axs[1])
    cbar2.set_label('Z Value')
    axs[1].set_xlabel('Index ID')
    axs[1].set_title('WT-MUT Matrix Plot')

    # Third scatter plot
    sc3 = axs[2].scatter(sub['Index1'], sub['Index2'], c=sub['Delta_Distance'], cmap='viridis')
    cbar3 = fig.colorbar(sc3, ax=axs[2])
    cbar3.set_label('Z Value')
    axs[2].set_xlabel('Index ID')
    axs[2].set_title('WT-MUT Matrix Plot')

    # Adjust layout for better spacing between plots
    plt.tight_layout()

    buffer = io.BytesIO()
    plt.savefig(buffer, format="png")
    buffer.seek(0)
    salt_plot_image = base64.b64encode(buffer.getvalue()).decode('utf8')
    buffer.close()
    plt.close()

    sub = sub[(sub['Distance_wt'] < 15) & (sub['Distance_mut'] < 15)]

    sub_rm=sub.query('ResidueName1 != ResidueName2')
    sub_rm['Res-pair']=sub_rm['ResidueName1'].astype(str)+'-'+sub_rm['ResidueName2'].astype(str)
    sub_rm = sub_rm.query('`Res-pair` != "ARG-LYS" and `Res-pair` != "ASP-GLU"')
    sub_rm = sub_rm.query('`Res-pair` != "LYS-ARG" and `Res-pair` != "GLU-ASP"')

    sub_rm.drop(['Res-pair'],axis=1,inplace=True)

    delta_distances = sub_rm['Delta_Distance']

    # Calculate the number of bins using the square root choice
    num_bins = int(np.ceil(np.sqrt(len(delta_distances))))

    # Plot the distribution of filtered distances
    plt.hist(delta_distances, bins=num_bins, edgecolor='black')
    plt.xlabel('Distance')
    plt.ylabel('Frequency')
    plt.title(f'Distribution of Filtered Distances')

    buffer = io.BytesIO()
    plt.savefig(buffer, format="png")
    buffer.seek(0)
    distribution_graph = base64.b64encode(buffer.getvalue()).decode('utf8')
    buffer.close()
    plt.close()

    
    sub_rm.to_csv('Subtract_Files/salt_bridge.csv',index=False)

    return salt_plot_image, distribution_graph