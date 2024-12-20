import { React, useState, useEffect } from 'react';
import { BrowserRouter as Router, Route, Routes, Link, useLocation } from 'react-router-dom';
import Home from './Home';
import Allosteric from './Allosteric'; 
import './App.css';
import axios from 'axios';
import("3dmol/build/3Dmol.js").then( ($3Dmol) => {
    console.log($3Dmol);
    //can do things with $3Dmol here
    });

import MyImage from "./assets/Subtract-intro.png";

function PDBCompare() {

    const [pdbFile1, setPdbFile1] = useState(null);
    const [pdbFile2, setPdbFile2] = useState(null);
    const [pdbField1, setPdbField1] = useState('');
    const [pdbField2, setPdbField2] = useState('');
    const [edgeFile, setEdgeFile] = useState(null);
    const [subtractionPlot, setsubtractionPlot] = useState(null);
    const [saltPlot, setSaltPlot] = useState(null);
    const [distributionPlot, setDistributionPlot] = useState(null);
    const [subtractDistributionPlot, setSubtractDistributionPlot] = useState(null);
    const [saltDistributionPlot, setSaltDistributionPlot] = useState(null);
    const [activePlotTypeTab, setActivePlotTypeTab] = useState(0);
    const [activePDBUploadTab, setActivePDBUploadTab] = useState(0);
    const [showRenderOptions, setShowRenderOptions] = useState(false);
    const [lowerBound, setLowerBound] = useState('');
    const [isLoading, setIsLoading] = useState(false);
    const [selectedChains, setSelectedChains] = useState({});
    const [chainRanges, setChainRanges] = useState({});
    const [edgesTable, setEdgesTable] = useState([]);
    const [uniqueID, setUniqueID] = useState('');

    const handleRerender = async () => {
        const formData = new FormData();
        formData.append('unique_id', uniqueID)
        formData.append('lower_bound', lowerBound);
        formData.append('selected_chains', JSON.stringify(selectedChains));
        formData.append('chain_ranges', JSON.stringify(chainRanges));
        console.log(uniqueID, " is unique id");
        
        // Only append edgeFile if it exists
        if (edgeFile) {
            formData.append('edge_file', edgeFile);
        }

        try {
            let response;
            if(activePlotTypeTab == 0) { // subtract edges
                response = await axios.post('/api/rerender', formData, {
                    headers: {
                        'Content-Type': 'multipart/form-data',
                    },
                });
            }
            else { // salt edges
                response = await axios.post('/api/rerender-salt', formData, {
                    headers: {
                        'Content-Type': 'multipart/form-data',
                    },
                });
            }
            const data = response.data;
            const parsedTable = JSON.parse(data.table);
            console.log('Parsed table data:', parsedTable);

            setEdgesTable(parsedTable);

            if (data.edges.length === 0) {
                alert("No edges found in the data. Consider lowering delta distance cutoff");
            }

            render3dmol(data);

        } catch (error) {
            console.error('Error:', error);
            alert('An error occurred while processing the PDB files.');
        }
    };

    const switchPlotTypeTab = (tabIndex) => {
        setActivePlotTypeTab(tabIndex);
        if(tabIndex == 0) {
            setDistributionPlot(subtractDistributionPlot);
        }
        else {
            setDistributionPlot(saltDistributionPlot);
        }
    };

    const switchPDBUploadTab = (tabIndex) => {
        setActivePDBUploadTab(tabIndex);
    };

    const handlePdbFile1Change = (event) => {
        setPdbFile1(event.target.files[0]);
        setPdbField1('');
        setPdbField2('');
    };

    const handlePdbFile2Change = (event) => {
        setPdbFile2(event.target.files[0]);
    };

    const handleEdgesFileChange = (event) => {
        setEdgeFile(event.target.files[0]);
    };

    const handleChainChange = (chain) => {
        setSelectedChains((prevChains) => ({
            ...prevChains,
            [chain]: !prevChains[chain],
        }));
    };

    const handleRangeInput = (chain, input) => {
        const rangesArray = input.split(',').map(range => range.trim());
        const parsed = rangesArray.map(range => {
            const [start, end] = range.split('-').map(Number);
            return { start, end };
        });

        setChainRanges(prev => ({ ...prev, [chain]: input }));

    };

    const handleSubmit = async () => {
        if ((!pdbFile1 || !pdbFile2) && (pdbField1 === '' || pdbField2 === '')) {
            alert('Please select both PDB files.');
            return;
        }

        const formData = new FormData();
        if(activePDBUploadTab === 0) {
            formData.append('pdb_file1', pdbFile1);
            formData.append('pdb_file2', pdbFile2);
        } else {
            try {
                const response1 = await axios.get(`https://files.rcsb.org/download/${pdbField1}.pdb`);
                const response2 = await axios.get(`https://files.rcsb.org/download/${pdbField2}.pdb`);
        
                formData.append('pdb_file1', new Blob([response1.data], { type: 'text/plain' }), `${pdbField1}.pdb`);
                formData.append('pdb_file2', new Blob([response2.data], { type: 'text/plain' }), `${pdbField2}.pdb`);
            } catch (error) {
                console.error('Error fetching PDB files:', error);
                alert('Could not fetch one or both PDB files. Please check the PDB IDs.');
            }
        }

        try {
            const chain_response = await axios.post('/api/extract_chains', formData, {
                headers: {
                    'Content-Type': 'multipart/form-data',
                },
            });

            const initialSelectedChains = {};
            console.log(chain_response.data, "Chain response from backend");
            const chain_data = chain_response.data;
            const chains = chain_data.chains;
            chains.forEach(chain => {
                initialSelectedChains[chain] = true;
            });

            setSelectedChains(initialSelectedChains);
            console.log(initialSelectedChains, " are the selected chains");
        } catch (error) {
            console.error('Error:', error);
            alert('Chain Format Issue');
        }

        setIsLoading(true);
        try {
            const response = await axios.post('/api/subtract', formData, {
                headers: {
                    'Content-Type': 'multipart/form-data',
                },
            });
            const data = response.data;
            const plots = response.data;
            setsubtractionPlot(plots.calculated_matrix_image);
            setSaltPlot(plots.salt_plot_image);
            setSubtractDistributionPlot(plots.subtract_distribution_graph);
            setSaltDistributionPlot(plots.salt_distribution_graph);
            setUniqueID(data.unique_id);

            render3dmol(data);
            setShowRenderOptions(true);

        } catch (error) {
            console.error('Error:', error);
            alert('An error occurred while processing the PDB files.');
        } finally {
            setIsLoading(false);
        }
    };

    const render3dmol = async (data) => {
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

        data.edges.forEach(edge => {
            viewer.addCylinder(
              {start: {x: edge.coords.start[0], y: edge.coords.start[1], z: edge.coords.start[2]},
              end: {x: edge.coords.end[0], y: edge.coords.end[1], z: edge.coords.end[2]},
              radius: 0.5,
              color: "blue",
              hoverable: true,
              opacity: 0.9,
              hover_callback: function(atom, viewer, event, container) {
                // Show the tooltip when hovering
                tooltip.style.display = 'block';
                tooltip.style.left = `${event.clientX}px`;  // Position tooltip near mouse cursor
                tooltip.style.top = `${event.clientY + window.scrollY}px`;

                // Set the tooltip content with edge label
                tooltip.innerHTML = `Edge Label: ${edge.label}`;
              },
              unhover_callback: function(atom, viewer, event, container) {
                // Hide the tooltip when not hovering
                tooltip.style.display = 'none';
              }
             });
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

        viewer.zoomTo();                                      
        viewer.render();                                     
        viewer.zoom(1.2, 1000);   
    }

    useEffect(() => {
        if (!distributionPlot && activePlotTypeTab == 0) {
            setDistributionPlot(subtractDistributionPlot);
        } else if (!distributionPlot && activePlotTypeTab == 1) {
            setDistributionPlot(saltDistributionPlot);
        }
        console.log(selectedChains, " are the selected chains");
    }, [subtractDistributionPlot, distributionPlot, chainRanges]);

    return (
    <div>
        <h1>PDB Pair Compare</h1>
        <img src={MyImage} alt="Residue-Residue Distance" style={{ maxWidth: "100%", height: "auto", marginTop: "20px" }} />
        <br></br>

        <div className="tab-navigation">
            <button onClick={() => switchPDBUploadTab(0)} className={activePDBUploadTab === 0 ? 'red-tab' : ''}>Upload PDB's</button>
            <button onClick={() => switchPDBUploadTab(1)} className={activePDBUploadTab === 1 ? 'red-tab' : ''}>Fetch PDB's</button>
        </div>

        <div className="tab-content">
            {activePDBUploadTab === 0 && (
                <div>
                    <input type="file" onChange={handlePdbFile1Change} />
                    <input type="file" onChange={handlePdbFile2Change} />
                    <button onClick={handleSubmit}>Submit</button>
                </div>
            )}
            {activePDBUploadTab === 1 && (
                <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'center' }}>
                    Enter 2 PDB's Which Will Be Fetched from Protein Data Bank Database:
                    <div>
                        <input
                            type="text" // Set type to number for double input
                            value={pdbField1}
                            onChange={(e) => setPdbField1(e.target.value)}
                            style={{ marginLeft: '10px', padding: '5px' }}
                            placeholder="PDB 1" // Optional placeholder
                        />
                        <input
                            type="text" // Set type to number for double input
                            value={pdbField2}
                            onChange={(e) => setPdbField2(e.target.value)}
                            style={{ marginLeft: '10px', padding: '5px' }}
                            placeholder="PDB 2" // Optional placeholder
                        />
                        <button onClick={handleSubmit}>Submit</button>
                    </div>
                </div>
            )}

        </div>

        <br></br>

        {isLoading && <p>Loading, please wait...</p>}

        <div className="tab-navigation">
            <button onClick={() => switchPlotTypeTab(0)} className={activePlotTypeTab === 0 ? 'active-tab' : ''}>Pairwise Distances</button>
            <button onClick={() => switchPlotTypeTab(1)} className={activePlotTypeTab === 1 ? 'active-tab' : ''}>Salt Bridge Distances</button>
        </div>

        <div className="tab-content">
            {activePlotTypeTab === 0 && subtractionPlot && (
                <img src={`data:image/png;base64,${subtractionPlot}`} alt="Subtracted Matrix" className="centered-image" />
            )}
            {activePlotTypeTab === 1 && saltPlot && (
                <img src={`data:image/png;base64,${saltPlot}`} alt="Salt Plot" className="centered-image" />
            )}
        </div>

        {showRenderOptions && (
            <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'flex-start' }}>
                <label style={{ marginBottom: '5px' }}>
                    To Render Enter Delta Distance Value for Cutoff, Use Distribution Graph For Help:
                </label>
                {/* <div style={{ display: 'flex', alignItems: 'center' }}>
                    Optional, Submit Custom Edges CSV Table To Render: <input type="file" onChange={handleEdgesFileChange} />
                </div> */}
                <div style={{ display: 'flex', alignItems: 'center' }}>
                    <h3>Select Chains and Enter Ranges:</h3>
                    {Object.keys(selectedChains).map(chain => (
                        <div key={chain} style={{ marginBottom: '20px' }}>
                            <label>
                                <input
                                    type="checkbox"
                                    checked={selectedChains[chain]}
                                    onChange={() => handleChainChange(chain)}
                                />
                                Chain {chain}
                            </label>
                            {selectedChains[chain] && ( // Show range input only if the chain is selected
                                <div>
                                    <input
                                        type="text"
                                        value={chainRanges[chain]}
                                        onChange={(e) => handleRangeInput(chain, e.target.value)}
                                        style={{ marginLeft: '10px', padding: '5px' }}
                                        placeholder={`(e.g., 44-50, 100-110)`}
                                    />
                                </div>
                            )}
                        </div>
                    ))}
                </div>
                <input
                    type="number" // Set type to number for double input
                    value={lowerBound}
                    onChange={(e) => setLowerBound(e.target.value)}
                    step="0.01" // Set step to allow decimal values
                    style={{ marginLeft: '10px', padding: '5px' }}
                    placeholder="Delta Distance" // Optional placeholder
                />
                <button onClick={handleRerender} style={{ marginLeft: '10px', padding: '5px' }}>
                    Re-Render
                </button>
            </div>
        )}
        <div style={{ display: 'grid', gridTemplateColumns: '60% 40%', width: '100%'  }}>
            <div id="viewport" className="mol-container"></div>
            {distributionPlot && (
                <img 
                    src={`data:image/png;base64,${distributionPlot}`} 
                    alt="Distribution Graph" 
                    className="centered-image" 
                    style={{maxWidth: '100%', height: 'auto' }}
                />
            )}
        </div>
        <table>
            <thead>
                <tr>
                    <th>Dist-Pdb1</th>
                    <th>Dist-Pdb2</th>
                    <th>∆-Dist(Å)</th>
                    <th>ResID1</th>
                    <th>ChainID1</th>
                    <th>ResID2</th>
                    <th>ChainID2</th>
                </tr>
            </thead>
            <tbody>
                {edgesTable.map((row, index) => (
                    <tr key={index}>
                        <td>{row.Distance_wt.toFixed(1)}</td>
                        <td>{row.Distance_mut.toFixed(1)}</td>
                        <td>{row.Delta_Distance.toFixed(1)}</td>
                        <td>{row.ResidueID1}</td>
                        <td>{row.ChainID1}</td>
                        <td>{row.ResidueID2}</td>
                        <td>{row.ChainID2}</td>
                    </tr>
                ))}
            </tbody>
        </table>

    </div>

    );
}

export default PDBCompare;