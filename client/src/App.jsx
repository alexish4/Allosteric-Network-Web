import { useState, useEffect } from 'react';
import './App.css';
import * as d3 from 'd3';
import html2pdf from 'html2pdf.js';
import axios from 'axios';
// import 'jquery';
// import 'nglview-js-widgets';


function App() {
  const [count, setCount] = useState(0);
  const [file, setSelectedFile] = useState(null);
  const [tab, setTab] = useState('TopPaths');

  //variables for handleSubmit
  const [loading, setLoading] = useState(false);
  const [sourceResid, setSourceResid] = useState('');
  const [sinkResid, setSinkResid] = useState('');
  const [kResid, setKResid] = useState('');
  const [average, setAverage] = useState('Yes');
  const [startingIndexValue, setStartingIndexValue] = useState('0');
  const [imgData, setImgData] = useState('');
  const [imgData2, setImgData2] = useState('');
  const [paths, setPaths] = useState([]);
  const [paths2, setPaths2] = useState([]);
  const [nglContent, setNglContent] = useState(null);

  window.displayTopPaths = (paths, paths2) => {
    setPaths(paths);
    setPaths2(paths2);
  }

  function highlightPathEdges(path) {
    // Convert the path's nodes into a set of edges
    const edges = [];
    for (let i = 0; i < path.nodes.length - 1; i++) {
        edges.push({
            source: path.nodes[i],
            target: path.nodes[i + 1]
        });
    }

    // Update the styles of the links to highlight the edges in the path
    d3.selectAll(".links line")
    .attr("stroke", d => {
        // Check if the current link is part of the highlighted path
        const isHighlighted = edges.some(edge =>
            (d.source.id === edge.source && d.target.id === edge.target) ||
            (d.source.id === edge.target && d.target.id === edge.source)
        );
        return isHighlighted ? "orange" : "#999";
    })
    .attr("stroke-width", d => {
        // Increase the stroke width for highlighted edges
        const isHighlighted = edges.some(edge =>
            (d.source.id === edge.source && d.target.id === edge.target) ||
            (d.source.id === edge.target && d.target.id === edge.source)
        );
        return isHighlighted ? 16 : 8;
    });
  }

  // Expose displayHistogram globally
  window.displayHistogram = (imgData, imgData2) => {
      setImgData(imgData);
      setImgData2(imgData2);
  };

  const handleFileChange = (event) => {
    setSelectedFile(event.target.files[0]);
  };

  const openTab = (tabName) => {
    setTab(tabName);

    if (tabName === 'NGLView') {
      console.log("fetch test");
      fetchNGLContent();
    }
  };

  const fetchNGLContent = async () => {
    const script = document.createElement('script');
    script.src = "../node_modules/ngl/dist/ngl.js";

    console.log("fetch test 2");

    script.onload = () => {
      console.log("fetch test 3");
      const stage = new NGL.Stage("viewport");
      try {
        axios.post('http://127.0.0.1:5000/nglview', {
          headers: {
              'Content-Type': 'multipart/form-data',
          },
        })
        .then(response => {
            const data = response.data;

            console.log("fetch test 4");

            stage.loadFile("rcsb://1CRN").then(function (component) {
              // Apply representations from Flask data
              data.representations.forEach(rep => {
                component.addRepresentation(rep.type, {
                  selection: rep.selection,
                  alpha: rep.alpha || 1.0
                });
              });
              component.autoView();
            });

            // Create a new shape for cylinders (edges)
            var shape = new NGL.Shape("edgeShape", { disableImpostor: true });

            // Iterate through the edge list and add cylinders for each edge
            data.edges.forEach(edge => {
              // Add a cylinder for each edge, connecting two residues (resID1 and resID2)
              shape.addCylinder(
                edge.coords.start,  // Start position of the cylinder
                edge.coords.end,    // End position of the cylinder
                edge.color,         // Cylinder color
                edge.radius         // Cylinder radius
              );
            });

            // Add the shape to the stage
            var shapeComp = stage.addComponentFromObject(shape);
            shapeComp.addRepresentation("buffer");

            // Adjust view to include the cylinders
            shapeComp.autoView();
        })
      } catch (error) {
        console.error('Error fetching NGL content:', error);
      }
    };

    document.body.appendChild(script);
  };

  // Load external libraries
  const loadScript = (src) => {
    const script = document.createElement('script');
    script.src = src;
    script.async = true;
    document.body.appendChild(script);
  };

  useEffect(() => {
    // Load the required external scripts
    loadScript("https://cdnjs.cloudflare.com/ajax/libs/html2pdf.js/0.9.2/html2pdf.bundle.min.js");

    // Assuming graph.js is now part of your src/ directory, import it like a module
    import('./graph.js').then((module) => {
        // graph.js can initialize here or include setup in a function in graph.js
        //module.initGraph(); // This assumes there's an initGraph function in graph.js
    });

    // Any additional JavaScript from the HTML should go inside useEffect for React
  }, []);


  return (
    <div className="App">
      <h1>Current-Flow-Allostery</h1>
      <form id="upload-form" encType="multipart/form-data">
          <div className="form-group">
              <label htmlFor="csv-file">Upload File:</label>
              <input type="file" id="csv-file" accept=".csv, .dat" />
              <button type="button" id="help-button">File Format Help</button>
          </div>
          <div className="radio-group">
              Label Starting Node to Start at 0 or 1?
              <label><input type="radio" name="option" defaultChecked value="0" /> 0</label>
              <label><input type="radio" name="option" value="1" /> 1</label>
          </div>
          <div className="radio-group2">
              Use Average Betweenness For Top Paths?
              <label><input type="radio" name="option2" defaultChecked value="Yes" /> Yes</label>
              <label><input type="radio" name="option2" value="No" /> No</label>
          </div>
          <label htmlFor="source-resid">Source Residue ID(s):</label>
          <input type="text" id="source-resid" placeholder="Enter source node(s), e.g., 1 or 1,2,3" />
          <br />
          <label htmlFor="sink-resid">Sink Residue ID(s):</label>
          <input type="text" id="sink-resid" placeholder="Enter sink node(s), e.g., 1 or 1,2,3" />
          <br />
          <label htmlFor="k-resid">Enter # of Top Paths (Optional)</label>
          <input type="text" id="k-resid" />
          <br />
          <button type="submit">Submit</button>
          <div id="loading-spinner" style={{ display: 'none' }}>
              <img src="../static/spinner.gif" alt="Loading... May take up to a few minutes" />
          </div>
      </form>
      <br />
      <div id="modal-overlay" className="modal-overlay"></div>

      <div id="modal" className="modal">
          <span id="modal-close" className="modal-close">&times;</span>
          <h3>File Upload Help</h3>
          <p>Select a CSV or DAT file to upload. This will be used for generating graphs and analyzing data. Correlation Values must be between 0 and 1!</p>
          <p>CSV Format:</p>
          <p>
              For CSVs, please use a 3-column structure in the following format: [resi], [resj], [Correlation Value of an Edge]. Use ',' as
              delimiter and include both: [0, 1, 0.5] and [1, 0, 0.5]. Don’t include column headers.
          </p>
          <p>DAT Format:</p>
          <p>Submit .dat files with a correlation matrix of the graph.</p>
      </div>

      <div id="label-container" style={{ textAlign: 'center' }}>
          <label id="source-node-label" style={{ display: 'none', color: 'green' }}>Source Node: </label>
          <br />
          <label id="sink-node-label" style={{ display: 'none', color: 'red' }}>Sink Node: </label>
          <br />
          <br />
          <label id="directions-label" style={{ display: 'none' }}>Choose 2 Resistors to Calculate Betweenness. Drag to Move</label>
      </div>

      <div id="graph-container" style={{ width: '100%', height: '600px', overflow: 'hidden', position: 'relative' }}>
          <svg id="graph-svg" width="100%" height="100%"></svg>
      </div>
      <button id="calculate-button" style={{ display: 'none' }}>Calculate</button>
      <button id="refresh-button" style={{ display: 'none' }}>Refresh</button>
      <button id="zoom-in-button">Zoom In</button>
      <button id="zoom-out-button">Zoom Out</button>
      <button id="reset-zoom-button">Reset Zoom</button>
      <button id="download-pdf" style={{ display: 'none' }}>Download to PDF</button>

      <pre id="output"></pre>

      <br />
      {/* Tab links */}
      <div className="tab">
        <button className="tablinks" onClick={() => openTab('TopPaths')}>Top Paths</button>
        <button className="tablinks" onClick={() => openTab('Histograms')}>Histograms</button>
        <button className="tablinks" onClick={() => openTab('CorrelationMatrix')}>Correlation Matrix</button>
        <button className="tablinks" onClick={() => openTab('RankedNodes')}>Ranked Nodes</button>
        <button className="tablinks" onClick={() => openTab('NGLView')}>NGL View</button>
      </div>

      {/* Tab content */}
      {tab === 'TopPaths' && (
        <div id="TopPaths" className="tabcontent">
          <h2>Top Shortest Paths</h2>
          <div style={{ display: "flex" }}>
            <div style={{ flex: 1 }}>
              <strong>Edge Length = -ln(betweenness)</strong>
              {paths.map((path, index) => (
                <div key={index} style={{ display: "flex", alignItems: "center", marginBottom: "10px" }}>
                  <div style={{ marginRight: "10px" }}>
                      Path {index + 1}: Total Path Length From Betweenness: {path.edge_length}, Path: {path.nodes.map(node => startingIndexValue === "1" ? node + 1 : node).join(" -> ")}
                  </div>
                  <button onClick={() => highlightPathEdges(path)}>Highlight</button>
                </div>
              ))}
            </div>
            <div style={{ flex: 1 }}>
              <strong>Edge Length = -ln(correlation)</strong>
              {paths2.map((path, index) => (
                <div key={index} style={{ display: "flex", alignItems: "center", marginBottom: "10px" }}>
                  <div style={{ marginRight: "10px" }}>
                      Path {index + 1}: Total Path Length From Correlation: {path.edge_length}, Path: {path.nodes.map(node => startingIndexValue === "1" ? node + 1 : node).join(" -> ")}
                  </div>
                  <button onClick={() => highlightPathEdges(path)}>Highlight</button>
                </div>
              ))}
            </div>
          </div>
        </div>
      )}

      {tab === 'Histograms' && (
        <div id="Histograms" className="tabcontent">
          <h2>Histograms</h2>
          <div id="histogram-container" style={{ textAlign: 'center', marginTop: '20px' }}>
              {imgData && <img src={`data:image/png;base64,${imgData}`} alt="Histogram 1" />}
          </div>
          <div id="histogram-container2" style={{ textAlign: 'center', marginTop: '20px' }}>
              {imgData2 && <img src={`data:image/png;base64,${imgData2}`} alt="Histogram 2" />}
          </div>
        </div>
      )}

      {tab === 'CorrelationMatrix' && (
        <div id="CorrelationMatrix" className="tabcontent">
          <svg id="correlation-svg"></svg>
        </div>
      )}

      {tab === 'RankedNodes' && (
        <div id="RankedNodes" className="tabcontent">
          <h2>Ranked Nodes</h2>
          <div id="ranked-nodes"></div>
        </div>
      )}
      {tab === 'NGLView' && (
        <div id="viewport" style={{ width: '100%', height: '400px' }}></div>
      )}
    </div>
  );
}

export default App;
