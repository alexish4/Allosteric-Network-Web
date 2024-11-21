import MDAnalysis as mda
import pandas as pd
from flask import request, jsonify
from Bio.PDB import PDBParser

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

def filter_by_residue_ids(df, common_residue_ids):
    """
    Filters the DataFrame to include only rows with Residue IDs that are in the common_residue_ids list.
    """
    # Filter the DataFrame for rows where the 'Residue ID' is in the list of common residue IDs
    filtered_df = df[df['Residue ID'].isin(common_residue_ids)]
    
    return filtered_df

def parse_ranges(ranges_str):
    ranges = []
    for range_part in ranges_str.split(','):
        range_part = range_part.strip()
        if '-' in range_part:
            try:
                start, end = map(int, range_part.split('-'))
                if start <= end:
                    ranges.append((start, end))
                else:
                    print(f"Invalid range: {range_part}")
            except ValueError:
                print(f"Invalid format: {range_part}")
    return ranges

def rerender_edgelist_from_mda_universe_and_residue_pairs(pubStrucUniverse, residue_pairs):
    edge_list = []

    for resID1, chainID1, resID2 , chainID2, distance in residue_pairs:
        residue1 = pubStrucUniverse.select_atoms(f"resid {resID1} and segid {chainID1} and name CA")
        residue2 = pubStrucUniverse.select_atoms(f"resid {resID2} and segid {chainID2} and name CA")

        empty_residue = False

        if not empty_residue:
            # Use MDAnalysis to calculate the center of mass for each residue
            crd1 = residue1.center_of_mass()
            crd2 = residue2.center_of_mass()

            resname1 = residue1.residues[0].resname
            resid1 = int(residue1.residues[0].resid)
            
            resname2 = residue2.residues[0].resname
            resid2 = int(residue2.residues[0].resid)
            
            # Create an edge label based on the residue names and IDs
            edgeLabel = f'∆ Distance: {distance} ({resID1}.{chainID1}-{resID2}.{chainID2})'

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

def create_residue_pairs_list(csv_file, filtered_chains, validated_ranges, lower_bound = 6.0):
    # Load the CSV file
    df = pd.read_csv(csv_file)
    
    # Determine cutoff
    filtered_df = df.loc[
        (df['Delta_Distance'] >= lower_bound) & 
        (df['ChainID1'].isin(filtered_chains)) & 
        (df['ChainID2'].isin(filtered_chains))
    ]

    if validated_ranges: # If user entered ranges
        # Create a mask for the ranges by iterating over each chain and its range(s)
        range_mask = False  # Start with an empty mask

        for chain, ranges in validated_ranges.items():
            for r in ranges:
                min_range, max_range = r
                # Combine range mask with OR condition to include all ranges
                range_mask |= (
                    ((filtered_df['ChainID1'] == chain) & 
                     (filtered_df['ResidueID1'] >= min_range) & 
                     (filtered_df['ResidueID1'] <= max_range)) |
                    ((filtered_df['ChainID2'] == chain) & 
                     (filtered_df['ResidueID2'] >= min_range) & 
                     (filtered_df['ResidueID2'] <= max_range))
                )

        # Include chains not in validated_ranges by using all residues for those chains
        remaining_chains = set(filtered_chains) - set(validated_ranges.keys())
        remaining_mask = (
            (filtered_df['ChainID1'].isin(remaining_chains)) |
            (filtered_df['ChainID2'].isin(remaining_chains))
        )
        filtered_df = filtered_df[range_mask | remaining_mask]


    residue_pairs = filtered_df[['ResidueID1', 'ChainID1', 'ResidueID2', 'ChainID2', 'Delta_Distance']].values.tolist()
    
    return residue_pairs, filtered_df

def extract_chains():
    try:
        pdb_file1 = request.files['pdb_file1']

        pdb_file1_path = 'Subtract_Files/pdb_file1.pdb'
        pdb_file1.save(pdb_file1_path)

        u = mda.Universe(pdb_file1_path)
            
        # Extract unique chain IDs
        chain_ids = u.atoms.segids  # Chain IDs are stored in 'segids'
        unique_chain_ids = sorted(set(chain_ids))  # Ensure they are unique and sorted

        print(unique_chain_ids)

        return jsonify({'chains': unique_chain_ids}), 200
    except Exception as e:
        # Handle errors gracefully and provide feedback
        return jsonify({'error': str(e)}), 500
