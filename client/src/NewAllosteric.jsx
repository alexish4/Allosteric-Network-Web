import { React, useState, useEffect } from 'react';
import './App.css';
import axios from 'axios';
import("3dmol/build/3Dmol.js").then( ($3Dmol) => {
    console.log($3Dmol);
    //can do things with $3Dmol here
    });


function NewAllosteric() {
    const [pdbFile, setPdbFile] = useState(null);
    const [datFile1, setDatFile1] = useState(null);
    const [datFile2, setDatFile2] = useState(null);
    const [flownessType, setFlownessType] = useState(0);
    const [activeSecondaryContentTab, setActiveSecondaryContentTab] = useState(0);
    const [sourceValues, setSourceValues] = useState('');
    const [sinkValues, setSinkValues] = useState('');
    const [numOfTopPaths, setNumOfTopPaths] = useState('');
    const [average, setAverage] = useState(0); 
    const [betweennessTopPaths1, setBetweennessTopPaths1] = useState([]);
    const [betweennessTopPaths2, setBetweennessTopPaths2] = useState([]);
    const [correlationTopPaths1, setCorrelationTopPaths1] = useState([]);
    const [correlationTopPaths2, setCorrelationTopPaths2] = useState([]);
    const [deltaValues, setDeltaValues] = useState(new Map());
    const [frequencyValues, setFrequencyValues] = useState(new Map());
    const [showResults, setShowResults] = useState(false);
    
    const [wtData, setWtData] = useState(null);
    const [mutData, setMutData] = useState(null);
    const [residueTable, setResidueTable] = useState([]);

    const handlePdbFileChange = (event) => {
        setPdbFile(event.target.files[0]);
    };

    const handleDatFile1Change = (event) => {
        setDatFile1(event.target.files[0]);
    };

    const handleDatFile2Change = (event) => {
        setDatFile2(event.target.files[0]);
    };

    const handleAverageChoice = (event) => {
        setAverage(parseInt(event.target.value));
    };

    const switchFlownessTypeTab = (tabIndex) => {
        setFlownessType(tabIndex);

        if (wtData !== null) {
            render3dmol(wtData, mutData, 0, tabIndex, residueTable, 0, deltaValues, frequencyValues);
        }
    };

    const switchSecondaryContentTab = (tabIndex) => {
        setActiveSecondaryContentTab(tabIndex);
    };

    const handleSubmit = async () => {
        const formData = new FormData();
        formData.append('pdb_file', pdbFile);
        formData.append('correlation_dat', datFile1);
        formData.append('source_values', sourceValues);
        formData.append('sink_values', sinkValues);
        formData.append('k', numOfTopPaths);
        formData.append('average', average);

        const formData2 = new FormData();
        formData2.append('pdb_file', pdbFile);
        formData2.append('correlation_dat', datFile2);
        formData2.append('source_values', sourceValues);
        formData2.append('sink_values', sinkValues);
        formData2.append('k', numOfTopPaths);
        formData2.append('average', average);

        try {
            const response = await axios.post('/api/allosteric', formData, {
                headers: {
                    'Content-Type' : 'multipart/form-data',
                },
            });

            const response2 = await axios.post('/api/allosteric', formData2, {
                headers: {
                    'Content-Type' : 'multipart/form-data',
                },
            });

            const wtData = response.data;
            const mutData = response2.data;
            const parsedTable = JSON.parse(wtData.table);
            setWtData(wtData);
            setMutData(mutData);
            setBetweennessTopPaths1(wtData.top_paths);
            setCorrelationTopPaths1(wtData.top_paths2);
            setBetweennessTopPaths2(mutData.top_paths);
            setCorrelationTopPaths2(mutData.top_paths2);
            setResidueTable(parsedTable);
            let [delta_map, frequencies_map] = calculateDeltaEdge(wtData.top_paths, mutData.top_paths);
            setDeltaValues(delta_map);
            setFrequencyValues(frequencies_map);
            render3dmol(wtData, mutData, 0, flownessType, parsedTable, -1, delta_map, frequencies_map); // by default don't highlight top path
            setShowResults(true);
        } catch (error) {
            console.error('Error:', error);
            alert('An error occurred while processing the files.');
        }        
    };

    const handleWtHighlight = (path, index) => {
        console.log('Highlight clicked for path:', path);
        console.log("Path index is: ", index);
        // Add your highlight logic here

        render3dmol(wtData, mutData, 0, flownessType, residueTable, index, deltaValues, frequencyValues);
    };

    const handleMutHighlight = (path, index) => {
        console.log('Highlight clicked for path:', path);
        // Add your highlight logic here

        render3dmol(wtData, mutData, 1, flownessType, residueTable, index, deltaValues, frequencyValues);
    };

    const calculateDeltaEdge = (paths1, paths2) => {
        const edgeFrequency1 = new Map();
        const edgeFrequency2 = new Map();
        const deltaMap = new Map();
    
        // Helper function to count edge frequencies
        const countEdgeFrequencies = (paths, edgeFrequency) => {
            paths.forEach((path) => {
                const nodes = path.nodes;
                for (let i = 0; i < nodes.length - 1; i++) {
                    const edge = `${nodes[i]}-${nodes[i + 1]}`; // Edge as "node1-node2"
                    edgeFrequency.set(edge, (edgeFrequency.get(edge) || 0) + 1);
                }
            });
        };
    
        // Count edge frequencies for paths1 and paths2
        countEdgeFrequencies(paths1, edgeFrequency1);
        countEdgeFrequencies(paths2, edgeFrequency2);
    
        // Calculate the delta for each unique edge
        const allEdges = new Set([
            ...edgeFrequency1.keys(),
            ...edgeFrequency2.keys()
        ]);
    
        let maxDelta = 0;
        let minDelta = 0;

        // First pass to compute deltas and find max absolute delta
        const rawDeltas = new Map();
        const frequenciesMap = new Map();
        allEdges.forEach((edge) => {
            const frequency1 = edgeFrequency1.get(edge) || 0;
            const frequency2 = edgeFrequency2.get(edge) || 0;
            const frequencies = `${frequency1} : ${frequency2}`;
            const delta = frequency1 - frequency2;
            rawDeltas.set(edge, delta);
            frequenciesMap.set(edge, frequencies);
            maxDelta = Math.max(maxDelta, delta);
            minDelta = Math.min(minDelta, delta);
        });        

        // Second pass to normalize deltas
        allEdges.forEach((edge) => {
            const delta = rawDeltas.get(edge);
            const deltaScale = delta < 0 ? minDelta : maxDelta;
            const normalizedDelta = deltaScale !== 0 ? (Math.abs(delta) / deltaScale) : 0; // Avoid division by zero
            deltaMap.set(edge, normalizedDelta);
        });

        return [deltaMap, frequenciesMap]; // Map with keys as "node1-node2" and values as normalized delta
    };    

    const render3dmol = async (wt_data, mut_data, graphIndex, flowType, parsedTable, top_path_index, deltaMap, frequencyMap) => {
        let universe = wt_data.pdb_content;
        let element = document.querySelector('#viewport');
        let config = { backgroundColor: 'white' };
        let viewer = $3Dmol.createViewer( element, config );
        viewer.addModel( universe, "pdb");  

        const tooltip = document.createElement('div');
        tooltip.style.position = 'absolute';
        tooltip.style.backgroundColor = '#fff';
        tooltip.style.border = '1px solid #ccc';
        tooltip.style.padding = '5px';
        tooltip.style.display = 'none';  // Hide by default
        document.body.appendChild(tooltip);

        console.log(graphIndex, " is graph index");

        // Determine edges based on flow type
        const edges = {
            wt: flowType === 0 ? wt_data.betweenness_edges : wt_data.correlation_edges,
            mut: flowType === 0 ? mut_data.betweenness_edges : mut_data.correlation_edges,
        };

        // Map to track highlighted edges
        const highlightedEdges = new Map();

        // Function to add cylinders
        const addCylinder = (edge, color, frequencies) => {
            viewer.addCylinder({
                start: { x: edge.coords.start[0], y: edge.coords.start[1], z: edge.coords.start[2] },
                end: { x: edge.coords.end[0], y: edge.coords.end[1], z: edge.coords.end[2] },
                radius: 0.5,
                color: color,
                hoverable: true,
                opacity: 1.0,
                hover_callback: function (atom, viewer, event, container) {
                    tooltip.style.display = "block";
                    tooltip.style.left = `${event.clientX}px`;
                    tooltip.style.top = `${event.clientY + window.scrollY}px`;
                    tooltip.innerHTML = `Edge Label: ${edge.label}, Wt-Mut Frequencies: ${frequencies}`;
                },
                unhover_callback: function () {
                    tooltip.style.display = "none";
                },
            });
        };

        // Function to convert RGB to HEX
        const rgbToHex = (r, g, b) => {
            const toHex = (value) => value.toString(16).padStart(2, '0');
            return `#${toHex(r)}${toHex(g)}${toHex(b)}`;
        };

        // Function to calculate color based on delta value
        const getColorFromDelta = (delta) => {
            let r = 0, g = 0, b = 0;

            if (delta < 0) {
                // Negative deltas are shades of blue
                b = 255
                const intensity = Math.min(255, Math.floor(255 * Math.abs(delta))); // Higher magnitude -> darker blue
                r = 255 - intensity
                g = 255 - intensity
            } else if (delta > 0) {
                // Positive deltas are shades of red
                r = 255
                const intensity = Math.min(255, Math.floor(255 * delta)); // Higher magnitude -> darker red
                b = 255 - intensity
                g = 255 - intensity
            } else {
                r = 255, g = 255, b = 255;
            }

            // Convert RGB to HEX and return the color
            return rgbToHex(r, g, b);
        };

        // Function to process and highlight edges
        const processEdges = (edges, isHighlight) => {
            edges.forEach((edge) => {
                const edgeKey = `${edge.node1_index}-${edge.node2_index}`;

                if(!deltaMap.has(edgeKey)) {
                    console.log("Edge ", edgeKey, " is not in delta map");
                }

                // Determine the edge color
                //let edgeColor;
                if (isHighlight && edge.path_index === top_path_index) {              
                    highlightedEdges.set(edgeKey, { ...edge, color: "orange" });
                    const frequency = frequencyMap.get(edgeKey) || "Empty";
                    addCylinder(edge, "orange", frequency);
                } 

                // Process the edge
                if (!highlightedEdges.has(edgeKey)) {
                    const delta = deltaMap.get(edgeKey) || 0; // Get delta or default to 0
                    const frequency = frequencyMap.get(edgeKey) || "Empty";
                    let edgeColor = getColorFromDelta(delta);
                    highlightedEdges.set(edgeKey, { ...edge, color: edgeColor });
                    addCylinder(edge, edgeColor, frequency);
                }
            });
        };

        // Highlight edges based on the selected graph
        processEdges(graphIndex === 0 ? edges.wt : edges.mut, true);

        // Render the remaining edges
        processEdges(edges.wt, false);
        processEdges(edges.mut, false);

        const model = viewer.getModel();
        let atoms = model.selectedAtoms({});
        let chains = new Set(atoms.map(atom => atom.chain));
        console.log("Chains detected in model:", chains);

        viewer.setHoverable({}, true,
            function (atom, viewer, event, container) {
              if (!atom.label) {
                  atom.label = viewer.addLabel(atom.resn + "." + atom.resi + "." + atom.chain, { position: atom, backgroundColor: 'mintcream', fontColor: 'black' });
              }
            },
            function (atom) {
              if (atom.label) {
                  viewer.removeLabel(atom.label);
                  delete atom.label;
              }
            }
          );

        const colors = ['#ffd3b6', '#ff8b94', '#dcedc1', '#ffffba', '#ffd4e5', '#eecbff', '#d4ffea', '#dbdcff', '#f7e7b4', '#ae5a41', '#68c4af', '#5a5255'];

        let colorIndex = 0;

        chains.forEach(chainID => {
            let color = colors[colorIndex % colors.length]; // Cycle through colors
            viewer.setStyle({ chain: chainID }, { cartoon: { color: color } });
            colorIndex++;
        });

        wt_data.source_values.forEach((source) => {
            const row = parsedTable.find((row) => row.NewIndex === source);

            if (row) {
                // Convert Chain ID (e.g., PROA to A, PROB to B, etc.)
                const chain = row["Chain ID"].replace("PRO", "");
        
                // Set the style dynamically using values from the parsedTable
                viewer.setStyle(
                    {
                        chain: chain, // Extracted and transformed chain
                        resi: row["Residue ID"], // Residue ID from the table
                        atom: row["Atom Name"] // Atom Name from the table
                    },
                    {
                        sphere: { radius: 2, color: 'red' } // Style for the selected atoms
                    }
                );
            } else {
                console.warn(`No matching row found in parsedTable for NewIndex: ${source}`);
            }
        });

        wt_data.sink_values.forEach((sink) => {
            const row = parsedTable.find((row) => row.NewIndex === sink);

            if (row) {
                // Convert Chain ID (e.g., PROA to A, PROB to B, etc.)
                const chain = row["Chain ID"].replace("PRO", "");
        
                // Set the style dynamically using values from the parsedTable
                viewer.setStyle(
                    {
                        chain: chain, // Extracted and transformed chain
                        resi: row["Residue ID"], // Residue ID from the table
                        atom: row["Atom Name"] // Atom Name from the table
                    },
                    {
                        sphere: { radius: 2, color: 'green' } // Style for the selected atoms
                    }
                );
            } else {
                console.warn(`No matching row found in parsedTable for NewIndex: ${sink}`);
            }
        });

        viewer.zoomTo();                                      
        viewer.render();                                     
        viewer.zoom(1.2, 1000);   
    };

    const ColorScale = () => {
        return (
            <svg width="500" height="60">
                {/* Gradient Definition */}
                <defs>
                    <linearGradient id="colorGradient" x1="0%" x2="100%" y1="0%" y2="0%">
                        <stop offset="0%" stopColor="blue" />
                        <stop offset="50%" stopColor="white" />
                        <stop offset="100%" stopColor="red" />
                    </linearGradient>
                </defs>
    
                {/* Gradient Bar */}
                <rect x="50" y="20" width="400" height="25" fill="url(#colorGradient)" stroke="black" />
    
                {/* Labels */}
                <text x="50" y="15" fontSize="14" fill="blue">More important in Mut</text>
                <text x="250" y="55" fontSize="14" fill="black" textAnchor="middle">Neutral</text>
                <text x="450" y="15" fontSize="14" fill="red" textAnchor="end">More important in Wt</text>
            </svg>
        );
    };    
    
    return (
        <div>
            <h1>Current-Flow-Allostery</h1>
            Please Submit PDB File: <input type="file" onChange={handlePdbFileChange} />
            <br></br>
            Please Submit DAT Files: <input type="file" onChange={handleDatFile1Change} /> <input type="file" onChange={handleDatFile2Change} />
            <br></br>

            Enter Source ID's
            <input
                type="text"
                value={sourceValues}
                onChange={(e) => setSourceValues(e.target.value)}
                style={{ marginLeft: '10px', padding: '5px' }}
                placeholder={`(e.g., 44-50, 100-110)`}
            />
            <br></br>
            Enter Sink ID's
            <input
                type="text"
                value={sinkValues}
                onChange={(e) => setSinkValues(e.target.value)}
                style={{ marginLeft: '10px', padding: '5px' }}
                placeholder={`(e.g., 44-50, 100-110)`}
            />
            <br></br>
            Enter Number Of Top Paths (Optional):
            <input
                type="number" // Set type to number for double input
                value={numOfTopPaths}
                onChange={(e) => setNumOfTopPaths(e.target.value)}
                style={{ marginLeft: '10px', padding: '5px' }}
                placeholder="Top K Paths" // Optional placeholder
            />
            <br></br>

            Use Average Betweenness?
            <input
                type="radio"
                id="YesOption"
                name="average"
                value="0"
                checked={average === 0}
                onChange={handleAverageChoice}
            />
            <label htmlFor="YesOption">Yes</label>

            <input
                type="radio"
                id="NoOption"
                name="average"
                value="1"
                checked={average === 1}
                onChange={handleAverageChoice}
            />
            <label htmlFor="NoOption">No</label>

            <button onClick={handleSubmit}>Submit</button>

            <div className="tab-navigation">
                <button onClick={() => switchFlownessTypeTab(0)} className={flownessType === 0 ? 'active-tab' : ''}>Betweenness</button>
                <button onClick={() => switchFlownessTypeTab(1)} className={flownessType === 1 ? 'active-tab' : ''}>Correlation</button>
            </div>
            
            {showResults && (
                <>
                    <div style={{ color: 'red', fontWeight: 'bold', marginTop: '5px' }}>
                        {sourceValues && `🔴 Source Values: ${sourceValues}`}
                    </div>

                    <div style={{ color: 'green', fontWeight: 'bold', marginTop: '5px' }}>
                        {sinkValues && `🟢 Sink Values: ${sinkValues}`}
                    </div>

                    <div style={{ textAlign: 'center', marginTop: '20px' }}>
                        <h3>Edges Colored By Frequency</h3>
                        <ColorScale />
                    </div>             
                </>
            )}

            <div id="viewport" className="mol-container"></div>

            <div className="tab-navigation">
                <button onClick={() => switchSecondaryContentTab(0)} className={activeSecondaryContentTab === 0 ? 'active-tab' : ''}>Top Paths</button>
                <button onClick={() => switchSecondaryContentTab(1)} className={activeSecondaryContentTab === 1 ? 'active-tab' : ''}>Residue LookUp</button>
            </div>

            <div className="tab-content">
                {activeSecondaryContentTab === 0 && (
                    <div style={{ display: 'flex', alignItems: 'center' }}>
                        <div>
                            <h3>Top Paths From Wt</h3>
                            <ol>
                                {betweennessTopPaths1.map((path, index) => (
                                    <li key={index}>
                                        <strong>Path Length:</strong> {path.edge_length} <br />
                                        <strong>Nodes:</strong> {path.nodes.join(' → ')}
                                        <button onClick={() => handleWtHighlight(path, index)}>Highlight</button>
                                    </li>
                                ))}
                            </ol>
                        </div>
                        <div>
                            <h3>Top Paths From Mut</h3>
                            <ol>
                                {betweennessTopPaths2.map((path, index) => (
                                    <li key={index}>
                                        <strong>Path Length:</strong> {path.edge_length} <br />
                                        <strong>Nodes:</strong> {path.nodes.join(' → ')}
                                        <button onClick={() => handleMutHighlight(path, index)}>Highlight</button>
                                    </li>
                                ))}
                            </ol>
                        </div>
                    </div>
                )}
                {activeSecondaryContentTab === 1 && (
                    <div>
                        <h3>Residue Table For Index LookUp</h3>
                        <table>
                            <thead>
                                <tr>
                                    <th>Index</th>
                                    <th>AtomName</th>
                                    <th>ResName</th>
                                    <th>ResID</th>
                                    <th>ChainID</th>
                                </tr>
                            </thead>
                            <tbody>
                                {residueTable.map((row, index) => (
                                    <tr key={index}>
                                        <td>{row.NewIndex}</td>
                                        <td>{row["Atom Name"]}</td>
                                        <td>{row["Residue Name"]}</td>
                                        <td>{row["Residue ID"]}</td>  
                                        <td>{row["Chain ID"]}</td>       
                                    </tr>
                                ))}
                            </tbody>
                        </table>
                    </div>
                )}

            </div>
            
        </div>

    );
}

export default NewAllosteric;