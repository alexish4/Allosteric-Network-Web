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

        try {
        const response = await axios.post('http://127.0.0.1:5000/subtract', formData, {
            headers: {
                'Content-Type': 'multipart/form-data',
            },
        });
        const plots = response.data;
        setPlot1(plots.calculated_matrix_image);
        setPlot2(plots.subtracted_distance_matrix_image);
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
    </div>
    );
}

export default Subtract;