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
    const [sourceValues, setSourceValues] = useState('');
    const [sinkValues, setSinkValues] = useState('');
    const [numOfTopPaths, setNumOfTopPaths] = useState('');
    const [average, setAverage] = useState(0); 
    const [betweennessTopPaths, setBetweennessTopPaths] = useState([]);
    const [correlationTopPaths, setCorrelationTopPaths] = useState([]);

    const handlePdbFileChange = (event) => {
        setPdbFile(event.target.files[0]);
    };

    const handleDatFileChange = (event) => {
        setDatFile(event.target.files[0]);
    }

    const handleAverageChoice = (event) => {
        setAverage(parseInt(event.target.value));
        console.log(average, " is average");
      };

    const switchGraphTypeTab = (tabIndex) => {
        setActiveGraphTypeTab(tabIndex);
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
            setBetweennessTopPaths(data.top_paths);
            setCorrelationTopPaths(data.top_paths2);
            console.log(data.top_paths);
            render3dmol(data);
        } catch (error) {
            console.error('Error:', error);
            alert('An error occurred while processing the files.');
        }
    };

    const handleHighlight = (path) => {
        console.log('Highlight clicked for path:', path);
        // Add your highlight logic here
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

            {/* For some reason Yes is counted as 1 and No is counted as 0 */}
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
            <div style={{ display: 'flex', alignItems: 'center' }}>
                <div>
                    <h3>Top Paths From Betweeness</h3>
                    <ol>
                        {betweennessTopPaths.map((path, index) => (
                            <li key={index}>
                                <strong>Edge Length:</strong> {path.edge_length} <br />
                                <strong>Nodes:</strong> {path.nodes.join(' → ')}
                                <button onClick={() => handleHighlight(path)}>Highlight</button>
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
                                <button onClick={() => handleHighlight(path)}>Highlight</button>
                            </li>
                        ))}
                    </ol>
                </div>
            </div>
        </div>

    );
}

export default NewAllosteric;