from flask import request, jsonify
import uuid
from flask import request, jsonify
import numpy as np
from scipy.sparse import coo_matrix
import copy
import networkx as nx
from itertools import islice
from collections import Counter
import logging
import MDAnalysis as mda
from PDBCompareMethods import pdb_to_dataframe

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

def process_dat_file(file):
    data = np.loadtxt(file)
    num_nodes = data.shape[0]
    rows = []
    cols = []
    correlations = []

    for i in range(num_nodes):
        for j in range(num_nodes):
            if i != j:  # Exclude self-loops
                mutual_info = data[i, j]
                if mutual_info > 0 and not np.isnan(mutual_info) and not np.isinf(mutual_info):  # Filter invalid weights:
                    # Add both (i, j) and (j, i) to ensure bidirectional edges
                    rows.extend([i, j])
                    cols.extend([j, i])
                    correlations.extend([mutual_info, mutual_info])

    return rows, cols, correlations

def parse_int_ranges(input_string):
  int_list = []
  for item in input_string.split(','):
    item = item.strip()
    if '-' in item:
      start, end = map(int, item.split('-'))
      int_list.extend(range(start, end + 1))
    else:
      int_list.append(int(item))
  return int_list

def process_graph_data():
    pdb_file = request.files['pdb_file']
    dat_file = request.files['correlation_dat']
    source_array = request.form['source_values']
    sink_array = request.form['sink_values']

    source_array = parse_int_ranges(source_array)
    sink_array = parse_int_ranges(sink_array)


    all = False #calculate average or all
    
    k = 51 #by default k is 10

    if request.form['k'] != "": #if k has input
        k = int(request.form['k'])
    if request.form['average'] == 0:
        all = True
    print(request.form['average'], "is average")
    print(all, "is all")
    
    if dat_file.filename.endswith('.dat'):
        rows, cols, correlations = process_dat_file(dat_file)
    else:
        return jsonify({'error': "Must Use Dat File!"}), 500
    
    # Generate unique filenames
    unique_id = uuid.uuid4().hex  # Generate a unique identifier
    pdb_file_path = f'Subtract_Files/{unique_id}_pdb_file1.pdb'
    pdb_file.save(pdb_file_path)
    
    pdb_df = pdb_to_dataframe(pdb_file_path)
    pdb_df = pdb_df.query('`Atom Name` == "CB" | (`Atom Name` == "CA" & `Residue Name` == "GLY")')
    pdb_df['NewIndex']=range(0,len(pdb_df)) # have indices match up with positions from correlation matrix
    print(pdb_df.head(), "is pdb head")
    print(len(pdb_df), "is length of df")

    adj_matrix = coo_matrix((correlations, (rows, cols)))

    # Create graph
    G = nx.from_scipy_sparse_array(adj_matrix)

    # label nodes to match their index in the adjacency matrix
    mapping = {i: i for i in range(adj_matrix.shape[0])}
    G = nx.relabel_nodes(G, mapping)

    tempAdjDense = adj_matrix.todense()
    
    # Convert adjacency matrix to Laplacian matrix
    tempLapDense = -copy.deepcopy(tempAdjDense)
    for ii, irow in enumerate(tempLapDense):
        tempLapDense[ii, ii] = -np.sum(irow)
    
    # Compute pseudoinverse of Laplacian matrix
    tempLinv = np.linalg.pinv(tempLapDense)

    #creating multiple graphs so you aren't just using average betweenness
    array_of_graphs = []
    if all:
        for so in source_array:
            for si in source_array:
                array_of_graphs.append(G.copy())
    print(len(array_of_graphs))

    #need the average graph because this is the graph we are drawing
    for u, v, data in G.edges(data=True):
        # Calculate based on the indices of the source (u) and target (v)
        betw = get_betw_value(u, v, tempLinv, tempAdjDense, source_array, sink_array)
        if betw is None:
            return jsonify({'error': "Empty Betweenness Score Found!"}), 500
        edge_length = -np.log(betw) #edge length is equal to -ln(|betw|) 
        edge_length2 = -np.log(data['weight'])

        data['betw'] = betw 
        data['edge_length'] = edge_length
        data['edge_length2'] = edge_length2

    if not all:
        top_paths, top_paths2, top_paths_lengths, top_paths2_lengths, most_important_nodes, most_important_nodes2 = generateTopPaths(G, k, source_array, sink_array)

    else:
        top_paths, top_paths2, top_paths_lengths, top_paths2_lengths, most_important_nodes, most_important_nodes2 = generateTopPaths2(array_of_graphs, k, tempLinv, tempAdjDense, source_array, sink_array)

    # Create the top_paths_data using path_lengths_edge_weights and top_paths_2
    print(top_paths, "is top paths")
    print(pdb_df.columns, "are columns")
    top_paths_data = [
        {'edge_length': top_paths_lengths[i], 'nodes': top_paths[i]}
        for i in range(len(top_paths))  
    ]
    top_paths_data2 = [
        {'edge_length': top_paths2_lengths[i], 'nodes': top_paths2[i]}
        for i in range(len(top_paths2))  
    ]

    ranked_nodes_data = [
        {'node': node, 'frequency': freq}
        for node, freq in most_important_nodes
    ]

    ranked_nodes_data2 = [
        {'node': node, 'frequency': freq}
        for node, freq in most_important_nodes2
    ]

    with open(pdb_file_path, 'r') as file:
        pdb_content = file.read()

    pdb_universe = mda.Universe(pdb_file_path)
    betweenness_edges = create_3d_edges(top_paths_data, pdb_df, G, pdb_universe)
    correlation_edges = create_3d_edges(top_paths_data2, pdb_df, G, pdb_universe)
    #edge_list = []

    graph_data = {
        'graph_data': nx.node_link_data(G),
        'top_paths': top_paths_data,
        'top_paths2': top_paths_data2,
        'ranked_nodes_data': ranked_nodes_data,
        'ranked_nodes_data2': ranked_nodes_data2,
        'pdb_content' : pdb_content,
        'betweenness_edges' : betweenness_edges,
        'correlation_edges' : correlation_edges,
        'table' : pdb_df.to_json(orient='records'),
        'unique_id' : unique_id
    }
    
    return graph_data

def create_3d_edges(top_paths_data, pdb_df, G, pdb_universe):
    edge_list = []

    for path in top_paths_data:
        nodes = path['nodes']
        for i in range(len(nodes) - 1):  # Pair adjacent nodes
            node1Index = nodes[i]
            node2Index = nodes[i+1]
            
            betweenness = G[node1Index][node2Index]["edge_length"]
            correlation = G[node1Index][node2Index]["edge_length2"]

            atom1 = pdb_df.query(f'NewIndex == {node1Index}').squeeze()  # Convert single-row DataFrame to Series
            atom2 = pdb_df.query(f'NewIndex == {node2Index}').squeeze()  # Convert single-row DataFrame to Series

            resID1 = atom1['Residue ID']
            resID2 = atom2['Residue ID']

            chainID1 = atom1['Chain ID']
            chainID2 = atom2['Chain ID']

            residue1 = pdb_universe.select_atoms(f"resid {resID1} and segid {chainID1} and name CA")
            residue2 = pdb_universe.select_atoms(f"resid {resID2} and segid {chainID2} and name CA")

            crd1 = residue1.center_of_mass()
            crd2 = residue2.center_of_mass()

            edgeLabel = f'Betweenness: {betweenness: .2f} Correlation: {correlation: .2f} ({resID1}.{chainID1}-{resID2}.{chainID2})'

            #converting to python types instead of numpy types so we can jsonify
            edge_data = {
                'label': edgeLabel,
                'coords': {
                    'start': [float(c) for c in crd1],  # Convert NumPy array to Python list of floats
                    'end': [float(c) for c in crd2]  # Convert NumPy array to Python list of floats
                }
            }
            
            edge_list.append(edge_data)
            
    return edge_list


    
def get_betw_value(u, v, tempLinv, tempAdjDense, source_array, sink_array):
    total_betweenness_score = 0

    try:
        for s in source_array:
            for t in sink_array:
                # Compute flow betweenness for the given edge
                v_source_sink_resist1 = tempLinv[u, s] - tempLinv[u, t]
                v_source_sink_resist2 = tempLinv[v, s] - tempLinv[v, t]
                b_resist1_resist2 = tempAdjDense[u, v] * (v_source_sink_resist1 - v_source_sink_resist2)
                total_betweenness_score += b_resist1_resist2

        #divide by number of combinations
        num_of_combinations = len(source_array) * len(sink_array)
        total_betweenness_score /= num_of_combinations

        betweenness_score = total_betweenness_score.item() # Convert to a standard Python type
        if betweenness_score < 0:
            betweenness_score *= -1

        return betweenness_score
    except IndexError as e:
        return None

def get_betw_value2(u, v, tempLinv, tempAdjDense, source, sink):
    # Compute flow betweenness for the given edge
    v_source_sink_resist1 = tempLinv[u, source] - tempLinv[u, sink]
    v_source_sink_resist2 = tempLinv[v, source] - tempLinv[v, sink]
    b_resist1_resist2 = tempAdjDense[u, v] * (v_source_sink_resist1 - v_source_sink_resist2)

    betweenness_score = b_resist1_resist2.item() # Convert to a standard Python type
    if betweenness_score < 0:
        betweenness_score *= -1

    return betweenness_score

def generateTopPaths(G, k, source_array, sink_array):
    top_paths = []
    top_paths_lengths = []
    top_paths2 = []
    top_paths2_lengths = []

    for so in source_array:
        for si in sink_array:
            # Find the top k optimal paths from source to sink
            paths = list(islice(nx.shortest_simple_paths(G, so, si, weight="edge_length"), k))
            paths2 = list(islice(nx.shortest_simple_paths(G, so, si, weight="edge_length2"), k))
            
            # Calculate path lengths for the first set of paths
            lengths = [sum(G[u][v]["edge_length"] for u, v in zip(path[:-1], path[1:])) for path in paths]
            lengths2 = [sum(G[u][v]["edge_length2"] for u, v in zip(path[:-1], path[1:])) for path in paths2]

            # Store the top paths and their lengths
            top_paths.extend(paths)
            top_paths_lengths.extend(lengths)
            top_paths2.extend(paths2)
            top_paths2_lengths.extend(lengths2)

    # Remove duplicates by converting paths to a tuple and using a set
    unique_paths_with_lengths = list({tuple(path): length for length, path in zip(top_paths_lengths, top_paths)}.items())
    unique_paths2_with_lengths = list({tuple(path): length for length, path in zip(top_paths2_lengths, top_paths2)}.items())
 
    # Sort the unique paths and lengths by length
    sorted_paths_with_lengths = sorted(unique_paths_with_lengths, key=lambda x: x[1])[:k]
    sorted_paths2_with_lengths = sorted(unique_paths2_with_lengths, key=lambda x: x[1])[:k]

    # Unpack the sorted pairs back into the arrays
    top_paths, top_paths_lengths = zip(*sorted_paths_with_lengths) if sorted_paths_with_lengths else ([], [])
    top_paths2, top_paths2_lengths = zip(*sorted_paths2_with_lengths) if sorted_paths2_with_lengths else ([], [])

    # Convert the results back to lists if needed
    top_paths = list(top_paths)
    top_paths_lengths = list(top_paths_lengths)
    top_paths2 = list(top_paths2)
    top_paths2_lengths = list(top_paths2_lengths)

    # Create an array of the most important nodes based on frequency
    all_nodes_in_paths = [node for path in top_paths for node in path]  # Flatten the list of top paths
    all_nodes_in_paths2 = [node for path in top_paths2 for node in path]  # Include nodes from the second set of paths
    node_frequencies = Counter(all_nodes_in_paths)  # Count the frequency of each node
    node_frequencies2 = Counter(all_nodes_in_paths2)

    # Create a list of (node, frequency_value) tuples sorted by frequency
    most_important_nodes = [(node, freq) for node, freq in node_frequencies.most_common()]
    most_important_nodes2 = [(node, freq) for node, freq in node_frequencies2.most_common()]

    print(most_important_nodes)
    print(most_important_nodes2)

    return top_paths, top_paths2, top_paths_lengths, top_paths2_lengths, most_important_nodes, most_important_nodes2

def generateTopPaths2(array_of_graphs, k, tempLinv, tempAdjDense, source_array, sink_array):
    array_index = 0
    for so in source_array:
        for si in sink_array:
            for u, v, data in array_of_graphs[array_index].edges(data=True):
                # Calculate based on the indices of the source (u) and target (v)
                epsilon = 1e-10 #prevent betw being so small that recorded as 0 which is bad for ln
                betw = get_betw_value2(u, v, tempLinv, tempAdjDense, so, si)
                edge_length = -np.log(max(betw, epsilon)) #edge length is equal to -ln(|betw|) 
                edge_length2 = -np.log(data['weight'])

                data['betw'] = betw 
                data['edge_length'] = edge_length
                data['edge_length2'] = edge_length2
            array_index += 1

    top_paths = []
    top_paths2 = []
    top_paths_lengths = []
    top_paths2_lengths = []

    array_index = 0
    for so in source_array:
        for si in sink_array:
            top_path, top_path2, length, length2 = miniGenerateTopPaths(array_of_graphs[array_index], k, so, si)
            top_paths.extend(top_path)
            top_paths_lengths.extend(length)
            top_paths2.extend(top_path2)
            top_paths2_lengths.extend(length2)
            array_index += 1
    #For betweenness
    unique_paths_with_lengths = list({tuple(path): length for length, path in zip(top_paths_lengths, top_paths)}.items())
    sorted_paths_with_lengths = sorted(unique_paths_with_lengths, key=lambda x: x[1])[:k]
    top_paths, top_paths_lengths = zip(*sorted_paths_with_lengths) if sorted_paths_with_lengths else ([], [])

    #for correlation
    unique_paths_with_lengths2 = list({tuple(path): length for length, path in zip(top_paths2_lengths, top_paths2)}.items())
    sorted_paths_with_lengths2 = sorted(unique_paths_with_lengths2, key=lambda x: x[1])[:k]
    top_paths2, top_paths2_lengths = zip(*sorted_paths_with_lengths2) if sorted_paths_with_lengths2 else ([], [])

    top_paths = list(top_paths)
    top_paths_lengths = list(top_paths_lengths)
    top_paths2 = list(top_paths2)
    top_paths2_lengths = list(top_paths2_lengths)

    # Create an array of the most important nodes based on frequency
    all_nodes_in_paths = [node for path in top_paths for node in path]  # Flatten the list of top paths
    all_nodes_in_paths2 = [node for path in top_paths2 for node in path]  # Include nodes from the second set of paths
    node_frequencies = Counter(all_nodes_in_paths)  # Count the frequency of each node
    node_frequencies2 = Counter(all_nodes_in_paths2)

    # Create a list of (node, frequency_value) tuples sorted by frequency
    most_important_nodes = [(node, freq) for node, freq in node_frequencies.most_common()]
    most_important_nodes2 = [(node, freq) for node, freq in node_frequencies2.most_common()]

    print(most_important_nodes)
    print(most_important_nodes2)

    return top_paths, top_paths2, top_paths_lengths, top_paths2_lengths, most_important_nodes, most_important_nodes2

def miniGenerateTopPaths(graph, k, so, si):
    top_paths = list(islice(nx.shortest_simple_paths(graph, so, si, weight="edge_length"), k))
    top_paths2 = list(islice(nx.shortest_simple_paths(graph, so, si, weight="edge_length2"), k))
    top_paths_lengths = [sum(graph[u][v]["edge_length"] for u, v in zip(path[:-1], path[1:])) for path in top_paths]
    top_paths2_lengths = [sum(graph[u][v]["edge_length2"] for u, v in zip(path[:-1], path[1:])) for path in top_paths2]

    return top_paths, top_paths2, top_paths_lengths, top_paths2_lengths