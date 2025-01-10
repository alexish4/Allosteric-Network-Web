import { React, useState, useEffect } from 'react';
import './App.css';
import axios from 'axios';
import("3dmol/build/3Dmol.js").then( ($3Dmol) => {
    console.log($3Dmol);
    //can do things with $3Dmol here
    });


function NewAllosteric() {
    const [pdbFile, setPdbFile] = useState(null);
    const [datFile, setDatFile] = useState(null);
    const [activeGraphTypeTab, setActiveGraphTypeTab] = useState(0);
    const [activeSecondaryContentTab, setActiveSecondaryContentTab] = useState(0);
    const [sourceValues, setSourceValues] = useState('');
    const [sinkValues, setSinkValues] = useState('');
    const [numOfTopPaths, setNumOfTopPaths] = useState('');
    const [average, setAverage] = useState(0); 
    const [betweennessTopPaths, setBetweennessTopPaths] = useState([]);
    const [correlationTopPaths, setCorrelationTopPaths] = useState([]);
    const [graphData, setGraphData] = useState(null);
    const [residueTable, setResidueTable] = useState([]);

    const handlePdbFileChange = (event) => {
        setPdbFile(event.target.files[0]);
    };

    const handleDatFileChange = (event) => {
        setDatFile(event.target.files[0]);
    }

    const handleAverageChoice = (event) => {
        setAverage(parseInt(event.target.value));
      };

    const switchGraphTypeTab = (tabIndex) => {
        setActiveGraphTypeTab(tabIndex);

        if (graphData !== null) {
            render3dmol(graphData, tabIndex, residueTable);
        }
    };

    const switchSecondaryContentTab = (tabIndex) => {
        setActiveSecondaryContentTab(tabIndex);
    };

    const handleSubmit = async () => {
        const formData = new FormData();
        formData.append('pdb_file', pdbFile);
        formData.append('correlation_dat', datFile);
        formData.append('source_values', sourceValues);
        formData.append('sink_values', sinkValues);
        formData.append('k', numOfTopPaths);
        formData.append('average', average);

        try {
            const response = await axios.post('/api/allosteric', formData, {
                headers: {
                    'Content-Type' : 'multipart/form-data',
                },
            });

            const data = response.data;
            const parsedTable = JSON.parse(data.table);
            setGraphData(data);
            setBetweennessTopPaths(data.top_paths);
            setCorrelationTopPaths(data.top_paths2);
            setResidueTable(parsedTable);
            render3dmol(data, activeGraphTypeTab, parsedTable, 0); // by default highlight top path
        } catch (error) {
            console.error('Error:', error);
            alert('An error occurred while processing the files.');
        }
    };

    const handleBetweennessHighlight = (path, index) => {
        console.log('Highlight clicked for path:', path);
        console.log("Path index is: ", index);
        // Add your highlight logic here

        render3dmol(graphData, 0, residueTable, index);
    };

    const handleCorrelationHighlight = (path, index) => {
        console.log('Highlight clicked for path:', path);
        // Add your highlight logic here

        render3dmol(graphData, 1, residueTable, index);
    };

    const render3dmol = async (data, graphIndex, parsedTable, top_path_index) => {
        let universe = data.pdb_content;
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

        let edges = graphIndex === 0 ? data.betweenness_edges : data.correlation_edges;
        const highlightedEdges = new Map();

        // render highlighted edges
        edges.forEach((edge) => {
            const edgeKey = `${edge.coords.start.join(",")}-${edge.coords.end.join(",")}`;
            const edgeColor = edge.path_index === top_path_index ? "blue" : "gray";
            if (edgeColor === "blue") {
                console.log("Updating edge to blue:", edgeKey);
                highlightedEdges.set(edgeKey, { ...edge, color: edgeColor });

                viewer.addCylinder({
                    start: { x: edge.coords.start[0], y: edge.coords.start[1], z: edge.coords.start[2] },
                    end: { x: edge.coords.end[0], y: edge.coords.end[1], z: edge.coords.end[2] },
                    radius: 0.5,
                    color: edgeColor,
                    hoverable: true,
                    opacity: 1.0,
                    hover_callback: function (atom, viewer, event, container) {
                        tooltip.style.display = "block";
                        tooltip.style.left = `${event.clientX}px`;
                        tooltip.style.top = `${event.clientY + window.scrollY}px`;
                        tooltip.innerHTML = `Edge Label: ${edge.label}`;
                    },
                    unhover_callback: function (atom, viewer, event, container) {
                        tooltip.style.display = "none";
                    },
                });
            }
        });

        // render the rest of the edges
        edges.forEach((edge) => {
            // Create a unique key for the edge based on start and end coordinates
            const edgeKey = `${edge.coords.start.join(",")}-${edge.coords.end.join(",")}`;
        
            if(!highlightedEdges.has(edgeKey)) {
                viewer.addCylinder({
                    start: { x: edge.coords.start[0], y: edge.coords.start[1], z: edge.coords.start[2] },
                    end: { x: edge.coords.end[0], y: edge.coords.end[1], z: edge.coords.end[2] },
                    radius: 0.5,
                    color: "gray",
                    hoverable: true,
                    opacity: 1.0,
                    hover_callback: function (atom, viewer, event, container) {
                        tooltip.style.display = "block";
                        tooltip.style.left = `${event.clientX}px`;
                        tooltip.style.top = `${event.clientY + window.scrollY}px`;
                        tooltip.innerHTML = `Edge Label: ${edge.label}`;
                    },
                    unhover_callback: function (atom, viewer, event, container) {
                        tooltip.style.display = "none";
                    },
                });
            }
        });

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

        const colors = ['orange', 'red', 'green', 'yellow', 'pink', 'purple', 'cyan', 'magenta', 'brown', 'silver', 'lime', 'black'];

        let colorIndex = 0;

        chains.forEach(chainID => {
            let color = colors[colorIndex % colors.length]; // Cycle through colors
            viewer.setStyle({ chain: chainID }, { cartoon: { color: color } });
            colorIndex++;
        });

        data.source_values.forEach((source) => {
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

        data.sink_values.forEach((sink) => {
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
    }

    return (
        <div>
            <h1>Current-Flow-Allostery</h1>
            Please Submit PDB File: <input type="file" onChange={handlePdbFileChange} />
            <br></br>
            Please Submit DAT File: <input type="file" onChange={handleDatFileChange} />
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
                <button onClick={() => switchGraphTypeTab(0)} className={activeGraphTypeTab === 0 ? 'active-tab' : ''}>Betweenness</button>
                <button onClick={() => switchGraphTypeTab(1)} className={activeGraphTypeTab === 1 ? 'active-tab' : ''}>Correlation</button>
            </div>


            <div id="viewport" className="mol-container"></div>

            <div className="tab-navigation">
                <button onClick={() => switchSecondaryContentTab(0)} className={activeSecondaryContentTab === 0 ? 'active-tab' : ''}>Top Paths</button>
                <button onClick={() => switchSecondaryContentTab(1)} className={activeSecondaryContentTab === 1 ? 'active-tab' : ''}>Residue LookUp</button>
            </div>

            <div className="tab-content">
                {activeSecondaryContentTab === 0 && (
                    <div style={{ display: 'flex', alignItems: 'center' }}>
                        <div>
                            <h3>Top Paths From Betweeness</h3>
                            <ol>
                                {betweennessTopPaths.map((path, index) => (
                                    <li key={index}>
                                        <strong>Edge Length:</strong> {path.edge_length} <br />
                                        <strong>Nodes:</strong> {path.nodes.join(' → ')}
                                        <button onClick={() => handleBetweennessHighlight(path, index)}>Highlight</button>
                                    </li>
                                ))}
                            </ol>
                        </div>
                        <div>
                            <h3>Top Paths From Correlation</h3>
                            <ol>
                                {correlationTopPaths.map((path, index) => (
                                    <li key={index}>
                                        <strong>Edge Length:</strong> {path.edge_length} <br />
                                        <strong>Nodes:</strong> {path.nodes.join(' → ')}
                                        <button onClick={() => handleCorrelationHighlight(path, index)}>Highlight</button>
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