import React, { useState } from 'react';
import './App.css';
import axios from 'axios';
let threeDMolLibrary = null;
import("3dmol/build/3Dmol.js").then((threeDMolModule) => {
    threeDMolLibrary = threeDMolModule.default ?? threeDMolModule;
    console.log(threeDMolLibrary);
    //can do things with $3Dmol here
    });


function NewAllosteric() {
    const [pdbFile1, setPdbFile1] = useState(null);
    const [pdbFile2, setPdbFile2] = useState(null);
    const [dcdFile1, setDCDFile1] = useState(null);
    const [dcdFile2, setDCDFile2] = useState(null);
    const [activeSystem, setActiveSystem] = useState(0);
    const [activePathMetric, setActivePathMetric] = useState(0);
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
    const [minDelta, setMinDelta] = useState(0);
    const [maxDelta, setMaxDelta] = useState(0);
    const [isLoading, setIsLoading] = useState(false);
    
    const [wtData, setWtData] = useState(null);
    const [mutData, setMutData] = useState(null);
    const [residueTable1, setResidueTable1] = useState([]);
    const [residueTable2, setResidueTable2] = useState([]);

    const FIGURE_URL = "Robust_Determination_of_Protein_Allosteric_Signaling_Pathways.png";

    const handlePdbFile1Change = (event) => {
        setPdbFile1(event.target.files[0]);
    };

    const handlePdbFile2Change = (event) => {
        setPdbFile2(event.target.files[0]);
    };

    const handleDCDFile1Change = (event) => {
        setDCDFile1(event.target.files[0]);
    };

    const handleDCDFile2Change = (event) => {
        setDCDFile2(event.target.files[0]);
    };

    const handleAverageChoice = (event) => {
        setAverage(parseInt(event.target.value));
    };

    const switchSystemTab = (systemIndex) => {
        setActiveSystem(systemIndex);

        if (wtData !== null) {
            const system1Paths = activePathMetric === 0 ? wtData.top_paths : wtData.top_paths2;
            const system2Paths = activePathMetric === 0 ? mutData.top_paths : mutData.top_paths2;
            const [deltaMap, frequenciesMap] = calculateDeltaEdge(system1Paths, system2Paths);
            const selectedResidueTable = systemIndex === 0 ? residueTable1 : residueTable2;

            setDeltaValues(deltaMap);
            setFrequencyValues(frequenciesMap);
            render3dmol(wtData, mutData, systemIndex, activePathMetric, selectedResidueTable, -1, deltaMap, frequenciesMap);
        }
    };

    const switchSecondaryContentTab = (tabIndex) => {
        setActiveSecondaryContentTab(tabIndex);
    };

    const unhighlight = () => {
        const selectedResidueTable = activeSystem === 0 ? residueTable1 : residueTable2;
        render3dmol(wtData, mutData, activeSystem, activePathMetric, selectedResidueTable, -1, deltaValues, frequencyValues);
    }

    const handleSubmit = async () => {
        const formData = new FormData();
        formData.append('pdb_file', pdbFile1);
        formData.append('render_pdb', pdbFile1);
        formData.append('trajectory', dcdFile1);
        formData.append('source_values', sourceValues);
        formData.append('sink_values', sinkValues);
        formData.append('k', numOfTopPaths);
        formData.append('average', average);

        const formData2 = new FormData();
        formData2.append('pdb_file', pdbFile2);
        formData2.append('render_pdb', pdbFile2);
        formData2.append('trajectory', dcdFile2);
        formData2.append('source_values', sourceValues);
        formData2.append('sink_values', sinkValues);
        formData2.append('k', numOfTopPaths);
        formData2.append('average', average);

        setIsLoading(true);

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
            const parsedTable1 = JSON.parse(wtData.table);
            const parsedTable2 = JSON.parse(mutData.table);
            setWtData(wtData);
            setMutData(mutData);
            setBetweennessTopPaths1(wtData.top_paths);
            setCorrelationTopPaths1(wtData.top_paths2);
            setBetweennessTopPaths2(mutData.top_paths);
            setCorrelationTopPaths2(mutData.top_paths2);
            setResidueTable1(parsedTable1);
            setResidueTable2(parsedTable2);
            const defaultPathMetric = 0;
            const system1Paths = wtData.top_paths;
            const system2Paths = mutData.top_paths;
            let [delta_map, frequencies_map] = calculateDeltaEdge(system1Paths, system2Paths);
            setActivePathMetric(defaultPathMetric);
            setDeltaValues(delta_map);
            setFrequencyValues(frequencies_map);
            const selectedResidueTable = activeSystem === 0 ? parsedTable1 : parsedTable2;
            render3dmol(wtData, mutData, activeSystem, defaultPathMetric, selectedResidueTable, -1, delta_map, frequencies_map); // by default don't highlight top path
            setShowResults(true);
        } catch (error) {
            console.error('Error:', error);
            alert('An error occurred while processing the files.');
        }        
        setIsLoading(false);
    };

    const handlePathHighlight = (path, index, systemIndex, pathMetric) => {
        console.log('Highlight clicked for path:', path);
        console.log("Path index is: ", index);
        const system1Paths = pathMetric === 0 ? wtData.top_paths : wtData.top_paths2;
        const system2Paths = pathMetric === 0 ? mutData.top_paths : mutData.top_paths2;
        const [deltaMap, frequenciesMap] = calculateDeltaEdge(system1Paths, system2Paths);
        const selectedResidueTable = systemIndex === 0 ? residueTable1 : residueTable2;

        setActivePathMetric(pathMetric);
        setDeltaValues(deltaMap);
        setFrequencyValues(frequenciesMap);
        render3dmol(wtData, mutData, systemIndex, pathMetric, selectedResidueTable, index, deltaMap, frequenciesMap);
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
    
        let max_delta = 0;
        let min_delta = 0;

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
            max_delta = Math.max(max_delta, delta);
            min_delta = Math.min(min_delta, delta);
        });        

        setMinDelta(min_delta);
        setMaxDelta(max_delta);

        // Second pass to normalize deltas
        allEdges.forEach((edge) => {
            const delta = rawDeltas.get(edge);
            const deltaScale = delta < 0 ? min_delta : max_delta;
            const normalizedDelta = deltaScale !== 0 ? (Math.abs(delta) / deltaScale) : 0; // Avoid division by zero
            deltaMap.set(edge, normalizedDelta);
        });

        return [deltaMap, frequenciesMap]; // Map with keys as "node1-node2" and values as normalized delta
    };    

    const render3dmol = async (wt_data, mut_data, graphIndex, flowType, parsedTable, top_path_index, deltaMap, frequencyMap) => {
        const selectedData = graphIndex === 0 ? wt_data : mut_data;
        let universe = selectedData.pdb_content;
        let element = document.querySelector('#viewport');
        let config = { backgroundColor: 'white' };
        if (!threeDMolLibrary) {
            throw new Error('3Dmol.js has not finished loading.');
        }

        let viewer = threeDMolLibrary.createViewer( element, config );
        viewer.addModel( universe, "pdb");  

        const tooltip = document.createElement('div');
        tooltip.style.position = 'absolute';
        tooltip.style.backgroundColor = '#fff';
        tooltip.style.border = '1px solid #ccc';
        tooltip.style.padding = '5px';
        tooltip.style.display = 'none';  // Hide by default
        tooltip.style.color = 'black';

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
                hover_callback: function (atom, viewer, event) {
                    tooltip.style.display = "block";
                    tooltip.style.left = `${event.clientX}px`;
                    tooltip.style.top = `${event.clientY + window.scrollY}px`;
                    tooltip.innerHTML = `Edge Label: ${edge.label}, System 1 / System 2 path frequencies: ${frequencies}`;
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
                if (!highlightedEdges.has(edgeKey) && !isHighlight) {
                    const delta = deltaMap.get(edgeKey) || 0; // Get delta or default to 0
                    const frequency = frequencyMap.get(edgeKey) || "Empty";
                    let edgeColor = getColorFromDelta(delta);
                    highlightedEdges.set(edgeKey, { ...edge, color: edgeColor });
                    addCylinder(edge, edgeColor, frequency);
                }
            });
        };

        const selectedEdges = graphIndex === 0 ? edges.wt : edges.mut;

        // Highlight the requested path and render the remaining edges for the selected system.
        processEdges(selectedEdges, true);
        processEdges(selectedEdges, false);

        const model = viewer.getModel();
        let atoms = model.selectedAtoms({});
        let chains = new Set(atoms.map(atom => atom.chain));
        console.log("Chains detected in model:", chains);

        viewer.setHoverable({}, true,
            function (atom, viewer) {
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

        selectedData.source_values.forEach((source) => {
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
                        sphere: { radius: 1.0, color: 'red' } // Style for the selected atoms
                    }
                );
            } else {
                console.warn(`No matching row found in parsedTable for NewIndex: ${source}`);
            }
        });

        selectedData.sink_values.forEach((sink) => {
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
                        sphere: { radius: 1.0, color: 'green' } // Style for the selected atoms
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
                <text x="50" y="15" fontSize="14" fill="blue">{minDelta}, System 2 enriched</text>
                <text x="250" y="15" fontSize="14" fill="black" textAnchor="middle">0</text>
                <text x="450" y="15" fontSize="14" fill="red" textAnchor="end">{maxDelta}, System 1 enriched</text>
            </svg>
        );
    };

    const selectedSystemName = activeSystem === 0 ? 'System 1' : 'System 2';
    const displayedBetweennessPaths = activeSystem === 0 ? betweennessTopPaths1 : betweennessTopPaths2;
    const displayedCorrelationPaths = activeSystem === 0 ? correlationTopPaths1 : correlationTopPaths2;
    const displayedResidueTable = activeSystem === 0 ? residueTable1 : residueTable2;

    const residueByNode = new Map(
        displayedResidueTable.map((row) => [Number(row.NewIndex), row])
    );

    const getNodeDistance = (firstNode, secondNode) => {
        const firstResidue = residueByNode.get(Number(firstNode));
        const secondResidue = residueByNode.get(Number(secondNode));

        if (!firstResidue || !secondResidue) {
            return null;
        }

        const firstCoordinates = [firstResidue.X, firstResidue.Y, firstResidue.Z].map(Number);
        const secondCoordinates = [secondResidue.X, secondResidue.Y, secondResidue.Z].map(Number);

        if (![...firstCoordinates, ...secondCoordinates].every(Number.isFinite)) {
            return null;
        }

        return Math.sqrt(
            firstCoordinates.reduce(
                (sum, coordinate, index) => sum + ((coordinate - secondCoordinates[index]) ** 2),
                0
            )
        );
    };

    const getTotalPathDistance = (path) => {
        let totalDistance = 0;

        for (let index = 0; index < path.nodes.length - 1; index++) {
            const distance = getNodeDistance(path.nodes[index], path.nodes[index + 1]);

            if (distance === null) {
                return null;
            }

            totalDistance += distance;
        }

        return totalDistance;
    };

    const renderPathNodes = (path) => {
        const totalPathDistance = getTotalPathDistance(path);

        return (
            <div>
                <div
                    style={{
                        display: 'flex',
                        flexWrap: 'wrap',
                        alignItems: 'center',
                        gap: '8px',
                        margin: '8px 0'
                    }}
                >
                    {path.nodes.map((node, position) => {
                        const residue = residueByNode.get(Number(node));
                        const nextNode = path.nodes[position + 1];
                        const distance = nextNode === undefined ? null : getNodeDistance(node, nextNode);

                        return (
                            <React.Fragment key={`${node}-${position}`}>
                                <span
                                    style={{
                                        display: 'inline-flex',
                                        flexDirection: 'column',
                                        padding: '7px 10px',
                                        border: '1px solid #b8c2cc',
                                        borderRadius: '7px',
                                        backgroundColor: '#f7fafc',
                                        lineHeight: 1.35
                                    }}
                                    title={residue ? `Node ${node}: ${residue['Residue Name']} ${residue['Residue ID']}, ${residue['Chain ID']}, ${residue['Atom Name']}` : `Node ${node}`}
                                >
                                    <strong>Node {node}</strong>
                                    {residue && (
                                        <span>
                                            {residue['Residue Name']} {residue['Residue ID']} · {residue['Chain ID']} · {residue['Atom Name']}
                                        </span>
                                    )}
                                </span>

                                {nextNode !== undefined && (
                                    <span
                                        title={`Straight-line distance between representative atoms in ${selectedSystemName}`}
                                        style={{ color: '#4a5568', whiteSpace: 'nowrap', fontWeight: 600 }}
                                    >
                                        — {distance === null ? 'distance unavailable' : `${distance.toFixed(2)} Å`} →
                                    </span>
                                )}
                            </React.Fragment>
                        );
                    })}
                </div>
                <div style={{ margin: '8px 0', color: '#2d3748' }}>
                    <strong>Total 3D structural distance:</strong>{' '}
                    {totalPathDistance === null ? 'unavailable' : `${totalPathDistance.toFixed(2)} Å`}
                </div>
            </div>
        );
    };

    const residueTableHeaderStyle = {
        border: '1px solid #cbd5e0',
        padding: '9px 12px',
        backgroundColor: '#edf2f7',
        textAlign: 'left',
        whiteSpace: 'nowrap'
    };

    const residueTableCellStyle = {
        border: '1px solid #e2e8f0',
        padding: '8px 12px',
        textAlign: 'left',
        whiteSpace: 'nowrap'
    };
    
    return (
        <div>
            <h1>Current-Flow-Allostery</h1>

            <div className="abstract-figure-container">
                <img
                    src={FIGURE_URL}
                    alt="Current-Flow-Allostery overview"
                    className="abstract-figure"
                />
            </div>
            <div>
                <strong>System 1:</strong> PDB file <input type="file" onChange={handlePdbFile1Change} />
                DCD trajectory <input type="file" onChange={handleDCDFile1Change} />
            </div>
            <div>
                <strong>System 2:</strong> PDB file <input type="file" onChange={handlePdbFile2Change} />
                DCD trajectory <input type="file" onChange={handleDCDFile2Change} />
            </div>

            Enter Source IDs
            <input
                type="text"
                value={sourceValues}
                onChange={(e) => setSourceValues(e.target.value)}
                style={{ marginLeft: '10px', padding: '5px' }}
                placeholder={`(e.g., 44-50, 100-110)`}
            />
            <br></br>
            Enter Sink IDs
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
                <button onClick={() => switchSystemTab(0)} className={activeSystem === 0 ? 'active-tab' : ''}>System 1</button>
                <button onClick={() => switchSystemTab(1)} className={activeSystem === 1 ? 'active-tab' : ''}>System 2</button>
            </div>
            
            {showResults && (
                <>
                    <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'center', marginTop: '10px' }}>
                        {/* Left: Source & Sink Values */}
                        <div style={{ flex: 1 }}>
                            <div style={{ color: 'red', fontWeight: 'bold', marginBottom: '5px' }}>
                                {sourceValues && `🔴 Source Values: ${sourceValues}`}
                            </div>

                            <div style={{ color: 'green', fontWeight: 'bold' }}>
                                {sinkValues && `🟢 Sink Values: ${sinkValues}`}
                            </div>
                        </div>

                        {/* Right: Color Scale */}
                        <div style={{ flex: 1, textAlign: 'center' }}>
                            <h3>Edges Colored by Δ Path Frequency (System 1 − System 2)</h3>
                            <ColorScale />
                        </div>
                    </div>
                    <div style={{ fontWeight: 'bold' }}>
                        Hover Over Edges for More Info
                    </div>
                    
                </>
            )}
            {isLoading && <p>Loading, please wait...</p>}

            <div id="viewport" className="mol-container"></div>
            <button onClick={() => unhighlight()}>Unhighlight</button>

            <div className="tab-navigation">
                <button onClick={() => switchSecondaryContentTab(0)} className={activeSecondaryContentTab === 0 ? 'active-tab' : ''}>Top Paths</button>
                <button onClick={() => switchSecondaryContentTab(1)} className={activeSecondaryContentTab === 1 ? 'active-tab' : ''}>Residue LookUp</button>
            </div>

            <div className="tab-content">
                {activeSecondaryContentTab === 0 && (
                    <div>
                        <h2>{selectedSystemName} Top Paths</h2>
                        <p>
                            System 1 is the first PDB/DCD pair and System 2 is the second PDB/DCD pair.
                            Both systems use the same source and sink residues. Betweenness-weighted and
                            correlation-weighted paths are shown together for the selected system. Distances
                            shown between neighboring nodes are straight-line distances between their
                            representative atoms in the selected system&apos;s PDB. The total 3D structural distance
                            is their sum; neither value is a betweenness or correlation path score.
                        </p>
                        <div
                            style={{
                                display: 'grid',
                                gridTemplateColumns: 'repeat(auto-fit, minmax(380px, 1fr))',
                                gap: '18px',
                                alignItems: 'start'
                            }}
                        >
                            <div
                                style={{
                                    maxHeight: '720px',
                                    overflowY: 'auto',
                                    padding: '12px',
                                    border: '1px solid #d8dee6',
                                    borderRadius: '8px'
                                }}
                            >
                                <h3>Betweenness-Weighted Paths</h3>
                                <ol style={{ paddingRight: '18px' }}>
                                    {displayedBetweennessPaths.map((path, index) => (
                                        <li key={index} style={{ marginBottom: '14px' }}>
                                            {renderPathNodes(path)}
                                            <button onClick={() => handlePathHighlight(path, index, activeSystem, 0)}>
                                                Highlight path
                                            </button>
                                        </li>
                                    ))}
                                </ol>
                            </div>
                            <div
                                style={{
                                    maxHeight: '720px',
                                    overflowY: 'auto',
                                    padding: '12px',
                                    border: '1px solid #d8dee6',
                                    borderRadius: '8px'
                                }}
                            >
                                <h3>Correlation-Weighted Paths</h3>
                                <ol style={{ paddingRight: '18px' }}>
                                    {displayedCorrelationPaths.map((path, index) => (
                                        <li key={index} style={{ marginBottom: '14px' }}>
                                            {renderPathNodes(path)}
                                            <button onClick={() => handlePathHighlight(path, index, activeSystem, 1)}>
                                                Highlight path
                                            </button>
                                        </li>
                                    ))}
                                </ol>
                            </div>
                        </div>
                    </div>
                )}
                {activeSecondaryContentTab === 1 && (
                    <div>
                        <h3>{selectedSystemName} Residue Lookup and Coordinates</h3>
                        <p>
                            Coordinates are in Ångströms and are used to calculate the straight-line distances shown between path nodes.
                        </p>
                        <div style={{ overflowX: 'auto' }}>
                        <table style={{ width: '100%', borderCollapse: 'collapse' }}>
                            <thead>
                                <tr>
                                    <th style={residueTableHeaderStyle}>Node</th>
                                    <th style={residueTableHeaderStyle}>Residue</th>
                                    <th style={residueTableHeaderStyle}>Chain</th>
                                    <th style={residueTableHeaderStyle}>Representative Atom</th>
                                    <th style={residueTableHeaderStyle}>X (Å)</th>
                                    <th style={residueTableHeaderStyle}>Y (Å)</th>
                                    <th style={residueTableHeaderStyle}>Z (Å)</th>
                                </tr>
                            </thead>
                            <tbody>
                                {displayedResidueTable.map((row, index) => (
                                    <tr key={index}>
                                        <td style={residueTableCellStyle}>{row.NewIndex}</td>
                                        <td style={residueTableCellStyle}>{row["Residue Name"]} {row["Residue ID"]}</td>
                                        <td style={residueTableCellStyle}>{row["Chain ID"]}</td>
                                        <td style={residueTableCellStyle}>{row["Atom Name"]}</td>
                                        <td style={residueTableCellStyle}>{Number(row.X).toFixed(3)}</td>
                                        <td style={residueTableCellStyle}>{Number(row.Y).toFixed(3)}</td>
                                        <td style={residueTableCellStyle}>{Number(row.Z).toFixed(3)}</td>
                                    </tr>
                                ))}
                            </tbody>
                        </table>
                        </div>
                    </div>
                )}

            </div>
            
        </div>

    );
}

export default NewAllosteric;
