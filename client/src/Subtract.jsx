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
    const [plot1, setPlot1] = useState(null);
    const [plot2, setPlot2] = useState(null);

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
        setPlot1(plots.calculated_matrix_image);
        setPlot2(plots.subtracted_distance_matrix_image);

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

        } catch (error) {
        console.error('Error:', error);
        alert('An error occurred while processing the PDB files.');
        }
    };

    return (
    <div>
        <h1>PDB File Uploader</h1>
        <h3>Please Submit 2 PDB Files:</h3>
        <input type="file" onChange={handlePdbFile1Change} />
        <input type="file" onChange={handlePdbFile2Change} />
        <button onClick={handleSubmit}>Submit</button>
        {plot1 && <img src={`data:image/png;base64,${plot1}`} alt="Calculated Matrix" className="centered-image" />}
        {plot2 && <img src={`data:image/png;base64,${plot2}`} alt="Subtracted Distance Matrix" className="centered-image" />}
        <div id="viewport" class="mol-container"></div>
    </div>
    );
}

export default Subtract;