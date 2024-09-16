import { useState, useEffect } from 'react';
import './App.css';
import * as d3 from 'd3';
import html2pdf from 'html2pdf.js';
import axios from 'axios';

function App() {
  const [count, setCount] = useState(0);
  const [selectedFile, setSelectedFile] = useState(null);
  const [tab, setTab] = useState('TopPaths');

  //variables for handleSubmit
  const [loading, setLoading] = useState(false);
  const [sourceResid, setSourceResid] = useState('');
  const [sinkResid, setSinkResid] = useState('');
  const [kResid, setKResid] = useState('');
  const [average, setAverage] = useState('');
  const [startingIndexValue, setStartingIndexValue] = useState('0');
  const [file, setFile] = useState(null);

  const handleFileChange = (event) => {
    setSelectedFile(event.target.files[0]);
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    setLoading(true); // Show loading spinner

    // Process source and sink inputs
    let sourceIdArray = sourceResid.split(',')
      .map(node => node.trim())
      .filter(node => node !== '');

    let sinkIdArray = sinkResid.split(',')
      .map(node => node.trim())
      .filter(node => node !== '');

    if (startingIndexValue === "1") {
      sourceIdArray = sourceIdArray.map(node => parseInt(node, 10) - 1); // Convert to integer
      sinkIdArray = sinkIdArray.map(node => parseInt(node, 10) - 1);     // Convert to integer
    } else {
      sourceIdArray = sourceIdArray.map(node => parseInt(node, 10)); // Convert to integer
      sinkIdArray = sinkIdArray.map(node => parseInt(node, 10));     // Convert to integer
    }

    const formData = new FormData();
    formData.append('file', file);
    formData.append('source', sourceIdArray);
    formData.append('sink', sinkIdArray);
    formData.append('k', kResid);
    formData.append('average', average);

    try {
      const response = await axios.post('http://your-server-url/upload', formData, {
        headers: {
          'Content-Type': 'multipart/form-data'
        }
      });

      const data = response.data;

      if (data.incorrect_input) {
        alert("Input Out of Bounds");
      } else {
        // Handle success: update source node label, display top paths, etc.
        updateSourceNodeLabel();
        displayTopPaths(data.top_paths, data.top_paths2); 
        displayRankedNodes(data.ranked_nodes_data, data.ranked_nodes_data2);
        drawGraph(data.graph_data);
        setupColorScaleAndEdges();
        drawColorScale();
        displayHistogram(data.histogram1, data.histogram2);
        drawCorrelationMatrix(data.graph_data);
      }

    } catch (error) {
      console.error('Error uploading file:', error);
    } finally {
      setLoading(false); // Hide loading spinner
    }
  };

  const openTab = (tabName) => {
    setTab(tabName);
  };

  useEffect(() => {
    if (tab === 'CorrelationMatrix') {
      // You can initialize D3.js code here
      const svg = d3.select('#correlation-svg');
      // D3.js code for drawing correlation matrix goes here
    }
  }, [tab]);

  return (
    <div className="App">
      <h1>Current-Flow-Allostery</h1>
      <form id="upload-form" onSubmit={handleSubmit}>
        <div className="form-group">
          <label htmlFor="csv-file">Upload File:</label>
          <input type="file" id="csv-file" accept=".csv, .dat" onChange={handleFileChange} />
          <button type="button" id="help-button">File Format Help</button>
        </div>
        <div className="radio-group">
          Label Starting Node to Start at 0 or 1?
          <label><input type="radio" name="option" value="0" defaultChecked /> 0</label>
          <label><input type="radio" name="option" value="1" /> 1</label>
        </div>
        <div className="radio-group2">
          Use Average Betweenness For Top Paths?
          <label><input type="radio" name="option2" value="Yes" defaultChecked /> Yes</label>
          <label><input type="radio" name="option2" value="No" /> No</label>
        </div>
        <label htmlFor="source-resid">Source Residue ID(s):</label>
        <input type="text" id="source-resid" placeholder="Enter source node(s), e.g., 1 or 1,2,3" />
        <br />
        <label htmlFor="sink-resid">Sink Residue ID(s):</label>
        <input type="text" id="sink-resid" placeholder="Enter sink node(s), e.g., 1 or 1,2,3" />
        <br />
        <label htmlFor="k-resid">Enter # of Top Paths(Optional)</label>
        <input type="text" id="k-resid" />
        <br />
        <button type="submit">Submit</button>
        <div id="loading-spinner" style={{ display: 'none' }}>
          <img src="../static/spinner.gif" alt="Loading..." />
        </div>
      </form>
      <br />
      {/* Tab links */}
      <div className="tab">
        <button className="tablinks" onClick={() => openTab('TopPaths')}>Top Paths</button>
        <button className="tablinks" onClick={() => openTab('Histograms')}>Histograms</button>
        <button className="tablinks" onClick={() => openTab('CorrelationMatrix')}>Correlation Matrix</button>
        <button className="tablinks" onClick={() => openTab('RankedNodes')}>Ranked Nodes</button>
      </div>

      {/* Tab content */}
      {tab === 'TopPaths' && (
        <div id="TopPaths" className="tabcontent">
          <h2>Top Shortest Paths</h2>
          <div id="top-paths"></div>
        </div>
      )}

      {tab === 'Histograms' && (
        <div id="Histograms" className="tabcontent">
          <h2>Histograms</h2>
          <div id="histogram-container" style={{ textAlign: 'center', marginTop: '20px' }}></div>
          <div id="histogram-container2" style={{ textAlign: 'center', marginTop: '20px' }}></div>
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

      <button id="zoom-in-button">Zoom In</button>
      <button id="zoom-out-button">Zoom Out</button>
      <button id="reset-zoom-button">Reset Zoom</button>
    </div>
  );
}

export default App;
