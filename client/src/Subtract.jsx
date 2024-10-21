// Import necessary components
import { React, useState } from 'react';
import { BrowserRouter as Router, Route, Routes, Link, useLocation } from 'react-router-dom';
import Home from './Home';
import Allosteric from './Allosteric'; 
import './App.css';
import axios from 'axios';
import("3dmol/build/3Dmol.js").then( ($3Dmol) => {
    console.log($3Dmol);
    //can do things with $3Dmol here
    });

function Subtract() {

    const [pdbFile1, setPdbFile1] = useState(null);
    const [pdbFile2, setPdbFile2] = useState(null);
    const [subtractionPlot, setsubtractionPlot] = useState(null);
    const [filteredPlot, setfilteredPlot] = useState(null);
    const [distributionPlot, setdistributionPlot] = useState(null);
    const [activeTab, setActiveTab] = useState(0);
    const [upperBound, setUpperBound] = useState('');
    const [lowerBound, setLowerBound] = useState('');

    // Handle button click
    const handleNewEnergyValue = async () => {
      console.log('Button clicked with input:', upperBound);
      if (!pdbFile1 || !pdbFile2) {
        alert('Please select both PDB files.');
        return;
        }

        const formData = new FormData();
        formData.append('pdb_file1', pdbFile1);
        formData.append('pdb_file2', pdbFile2);
        formData.append('upper_bound', upperBound);
        formData.append('lower_bound', lowerBound);

        try {
            const response = await axios.post('http://127.0.0.1:5000/rerender', formData, {
                headers: {
                    'Content-Type': 'multipart/form-data',
                },
            });
            const data = response.data;
            render3dmol(data);

        } catch (error) {
            console.error('Error:', error);
            alert('An error occurred while processing the PDB files.');
        }
    };

    // Function to handle tab switching
    const switchTab = (tabIndex) => {
        setActiveTab(tabIndex);
    };

    const handlePdbFile1Change = (event) => {
        setPdbFile1(event.target.files[0]);
    };

    const handlePdbFile2Change = (event) => {
        setPdbFile2(event.target.files[0]);
    };

    const handleSubmit = async () => {
        if (!pdbFile1 || !pdbFile2) {
        alert('Please select both PDB files.');
        return;
        }

        const formData = new FormData();
        formData.append('pdb_file1', pdbFile1);
        formData.append('pdb_file2', pdbFile2);

        console.log("Test")
        try {
            const response = await axios.post('http://127.0.0.1:5000/subtract', formData, {
                headers: {
                    'Content-Type': 'multipart/form-data',
                },
            });
            const data = response.data;
            const plots = response.data;
            setsubtractionPlot(plots.calculated_matrix_image);
            setfilteredPlot(plots.subtracted_distance_matrix_image);
            setdistributionPlot(plots.distribution_graph);

            render3dmol(data);

        } catch (error) {
            console.error('Error:', error);
            alert('An error occurred while processing the PDB files.');
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
                tooltip.style.top = `${event.clientY}px`;

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

        // Get all atoms in the model
        const atoms = model.selectedAtoms({});

        // Filter atoms by residue number and chain ID
        const residueAtoms = atoms.filter(atom => atom.resi === 543 && atom.chain === "A");

        // Print the coordinates of the filtered atoms
        residueAtoms.forEach(atom => {
            console.log(`Atom: ${atom.atom}, Chain: ${atom.chain}, Residue: ${atom.resi}, x: ${atom.x}, y: ${atom.y}, z: ${atom.z}`);
        });

        viewer.setStyle({chain:'A'}, {cartoon:{color:'orange'}});
        viewer.setStyle({chain:'B'}, {cartoon:{color:'red'}});
        viewer.setStyle({chain:'C'}, {cartoon:{color:'green'}});
        viewer.setStyle({chain:'D'}, {cartoon:{color:'yellow'}});
        viewer.zoomTo();                                      /* set camera */
        viewer.render();                                      /* render scene */
        viewer.zoom(1.2, 1000);   
    }

    return (
    <div>
        <h1>PDB File Uploader</h1>
        <h3>Please Submit 2 PDB Files:</h3>
        <input type="file" onChange={handlePdbFile1Change} />
        <input type="file" onChange={handlePdbFile2Change} />
        <button onClick={handleSubmit}>Submit</button>

        <div className="tab-navigation">
            <button onClick={() => switchTab(0)} className={activeTab === 0 ? 'active-tab' : ''}>Subtract Plot</button>
            <button onClick={() => switchTab(1)} className={activeTab === 1 ? 'active-tab' : ''}>Filtered Matrix</button>
            <button onClick={() => switchTab(2)} className={activeTab === 2 ? 'active-tab' : ''}>Distribution Graph</button>
        </div>

        <div className="tab-content">
            {activeTab === 0 && subtractionPlot && (
                <img src={`data:image/png;base64,${subtractionPlot}`} alt="Subtracted Matrix" className="centered-image" />
            )}
            {activeTab === 1 && filteredPlot && (
                <img src={`data:image/png;base64,${filteredPlot}`} alt="Filtered Matrix" className="centered-image" />
            )}
            {activeTab === 2 && distributionPlot && (
                <img src={`data:image/png;base64,${distributionPlot}`} alt="Distribution Graph" className="centered-image" />
            )}
        </div>

        <div style={{ display: 'flex', alignItems: 'center' }}>
            <div id="viewport" className="mol-container"></div>
            <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'flex-start' }}>
                <label style={{ marginBottom: '5px' }}>
                    To Re-render Enter New Lower and Upper Bounds For Energy Cutoff:
                </label>
                <input
                    type="number" // Set type to number for double input
                    value={lowerBound}
                    onChange={(e) => setLowerBound(e.target.value)}
                    step="0.01" // Set step to allow decimal values
                    style={{ marginLeft: '10px', padding: '5px' }}
                    placeholder="Lower Bound" // Optional placeholder
                />
                <input
                    type="number" // Set type to number for double input
                    value={upperBound}
                    onChange={(e) => setUpperBound(e.target.value)}
                    step="0.01" // Set step to allow decimal values
                    style={{ marginLeft: '10px', padding: '5px' }}
                    placeholder="Upper Bound" // Optional placeholder
                />
                <button onClick={handleNewEnergyValue} style={{ marginLeft: '10px', padding: '5px' }}>
                    Re-Render
                </button>
            </div>
        </div>
    </div>
    );
}

export default Subtract;