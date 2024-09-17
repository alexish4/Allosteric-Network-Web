import { useState, useEffect } from 'react';
import './App.css';
import * as d3 from 'd3';
import html2pdf from 'html2pdf.js';
import axios from 'axios';

function App() {
  const [count, setCount] = useState(0);
  const [file, setSelectedFile] = useState(null);
  const [tab, setTab] = useState('TopPaths');

  //variables for handleSubmit
  const [loading, setLoading] = useState(false);
  const [sourceResid, setSourceResid] = useState('');
  const [sinkResid, setSinkResid] = useState('');
  const [kResid, setKResid] = useState('');
  const [average, setAverage] = useState('');
  const [startingIndexValue, setStartingIndexValue] = useState('0');
  //const [file, setFile] = useState(null);

  let resistor_1 = null;
  let resistor_2 = null;

  let source_id_array = [];
  let sink_id_array = [];

  let large_bet = 0;

  const handleFileChange = (event) => {
    setSelectedFile(event.target.files[0]);
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    setLoading(true); // Show loading spinner

    console.log("Test");

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

    console.log(sourceIdArray);
    console.log(sinkIdArray);

    console.log(file);

    const formData = new FormData();
    formData.append('file', file);
    formData.append('source', sourceIdArray);
    formData.append('sink', sinkIdArray);
    formData.append('k', kResid);
    formData.append('average', average);

    try {
      const response = await axios.post('http://127.0.0.1:5000/upload', formData, {
        headers: {
          'Content-Type': 'multipart/form-data'
        }
      });

      const data = response.data;
      console.log("Test 2");
      console.log(data.largest_betweenness);

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
      <form onSubmit={handleSubmit} id="upload-form">
        Upload File:&nbsp;
        <input type="file" id="csv-file" accept=".csv, .dat" onChange={handleFileChange} required /><br></br>
        Source Residue ID/s&nbsp;
        <input type="text" id="source-resid" value={sourceResid} onChange={(e) => setSourceResid(e.target.value)} placeholder="Source Resid" required /><br></br>  
        Sink Residue ID/s&nbsp;
        <input type="text" id="sink-resid" value={sinkResid} onChange={(e) => setSinkResid(e.target.value)} placeholder="Sink Resid" required /><br></br>
        Enter # of Top Paths(Optional)&nbsp;
        <input type="text" id="k-resid" value={kResid} onChange={(e) => setKResid(e.target.value)} placeholder="K Resid" /><br></br>
        
        <div>
          <label>
            Label Starting Node to Start at 0 or 1? 
            <input type="radio" name="option" value="0" onChange={(e) => setStartingIndexValue(e.target.value)} checked={startingIndexValue === '0'} />
            0
          </label>
          <label>
            <input type="radio" name="option" value="1" onChange={(e) => setStartingIndexValue(e.target.value)} checked={startingIndexValue === '1'} />
            1
          </label>
        </div>

        <div>
          <label>
            Use Average Betweenness For Top Paths?
            <input type="radio" name="option2" value="Yes" onChange={(e) => setAverage(e.target.value)} checked={average === 'average1'} />
            Yes
          </label>
          <label>
            <input type="radio" name="option2" value="No" onChange={(e) => setAverage(e.target.value)} checked={average === 'average2'} />
            No
          </label>
        </div>

        <button type="submit">Upload</button>
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
