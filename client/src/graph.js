import axios from 'axios';
import * as d3 from 'd3';

let resistor_1 = null;
let resistor_2 = null;

let source_id_array = [];
let sink_id_array = [];

let large_bet = 0;
let startingIndexValue = 1;

const svg = d3.select("#graph-svg");
const width = parseInt(svg.style("width"));
const height = parseInt(svg.style("height"));

// Group for zoomable content
const zoomGroup = svg.append("g");

// Zoom and pan functionality
const zoom = d3.zoom()
    .scaleExtent([0.1, 10])
    .on("zoom", event => {
        zoomGroup.attr("transform", event.transform);
    });

svg.call(zoom);

// Zoom control buttons
d3.select("#zoom-in-button").on("click", () => zoom.scaleBy(svg.transition().duration(500), 1.2));
d3.select("#zoom-out-button").on("click", () => zoom.scaleBy(svg.transition().duration(500), 0.8));
d3.select("#reset-zoom-button").on("click", () => svg.transition().duration(500).call(zoom.transform, d3.zoomIdentity));

document.getElementById('help-button').addEventListener('click', function() {
    document.getElementById('modal').style.display = 'block';
    document.getElementById('modal-overlay').style.display = 'block';
});

document.getElementById('modal-close').addEventListener('click', function() {
    document.getElementById('modal').style.display = 'none';
    document.getElementById('modal-overlay').style.display = 'none';
});

document.getElementById('modal-overlay').addEventListener('click', function() {
    document.getElementById('modal').style.display = 'none';
    document.getElementById('modal-overlay').style.display = 'none';
});

document.getElementById('upload-form').onsubmit = function(e) {
    e.preventDefault();
    
    // Show the loading spinner
    document.getElementById('loading-spinner').style.display = 'block';

    const fileInput = document.getElementById('csv-file');
    const sourceResid = document.getElementById('source-resid').value;
    const sinkResid = document.getElementById('sink-resid').value;
    const k = document.getElementById('k-resid').value;

    //resetting resistor values
    resistor_1 = null
    resistor_2 = null

    // Get the selected radio button value
    startingIndexValue = document.querySelector('input[name="option"]:checked').value;
    if (startingIndexValue === "1") {
        console.log("Selected Starting Index:", startingIndexValue);
    }

    let average = document.querySelector('input[name="option2"]:checked').value;

    let sourceInput = document.getElementById('source-resid').value;
    let sinkInput = document.getElementById('sink-resid').value;

    // Process source input
    source_id_array = sourceInput
        .split(',')
        .map(node => node.trim())
        .filter(node => node !== '');

    if (source_id_array.length > 1) {
        document.getElementById("source-node-label").innerText = "Source Nodes: " + source_id_array.join(", ");
    } else{
        document.getElementById("source-node-label").innerText = "Source Node: " + source_id_array[0];
    }

    // Process sink input
    sink_id_array = sinkInput
        .split(',')
        .map(node => node.trim())
        .filter(node => node !== '');
        
    if (sink_id_array.length > 1) {
        document.getElementById("sink-node-label").innerText = "Sink Nodes: " + sink_id_array.join(", ");
    } else  {
        document.getElementById("sink-node-label").innerText = "Sink Node: " + sink_id_array[0];
    } 

    if (startingIndexValue === "1") {
        source_id_array = source_id_array.map(node => parseInt(node, 10) - 1); // Convert each to an integer
        sink_id_array = sink_id_array.map(node => parseInt(node, 10) - 1); // Convert each to an integer
    } else {
        source_id_array = source_id_array.map(node => parseInt(node, 10)); // Convert each to an integer
        sink_id_array = sink_id_array.map(node => parseInt(node, 10)); // Convert each to an integer
    }
    // Print arrays to console
    console.log("Source Nodes Array:", source_id_array);
    console.log("Sink Nodes Array:", sink_id_array);
    
    const formData = new FormData();
    formData.append('file', fileInput.files[0]);
    formData.append('source', source_id_array);
    formData.append('sink', sink_id_array);
    formData.append('k', k);
    formData.append('average', average);
    // Use axios instead of fetch
    axios.post('/api/upload', formData, {
        headers: {
            'Content-Type': 'multipart/form-data',
        },
    })
    .then(response => {
        const data = response.data;

        if (data.incorrect_input === true) {
            alert("Input Out of Bounds");
        } else {
            large_bet = data.largest_betweenness;
            updateSourceNodeLabel();
            window.displayTopPaths(data.top_paths, data.top_paths2); // Display the top paths
            displayRankedNodes(data.ranked_nodes_data, data.ranked_nodes_data2); // Display the ranked nodes
            drawGraph(data.graph_data); // Draw the graph
            setupColorScaleAndEdges();
            drawColorScale();
            window.displayHistogram(data.histogram1, data.histogram2); // Display histograms
            drawCorrelationMatrix(data.graph_data); // Draw the correlation matrix
        }
    })
    .catch(error => {
        console.error('Error uploading file:', error);
    })
    .finally(() => {
        // Hide the loading spinner after the process is complete
        document.getElementById('loading-spinner').style.display = 'none';
    });
};

function drawGraph(graph) {
    zoomGroup.selectAll("*").remove(); // Clear previous graph
    
    // Position the source and sink nodes
    const sourceNode = graph.nodes.find(d => d.id === source_id_array[0]);
    const sinkNode = graph.nodes.find(d => d.id === sink_id_array[0]);

    if (sourceNode) sourceNode.fx = 200; // Fixed x position for source node
    if (sinkNode) sinkNode.fx = width - 200; // Fixed x position for sink node

    const simulation = d3.forceSimulation(graph.nodes)
        .force("link", d3.forceLink(graph.links).id(d => d.id).distance(50).strength(1))
        .force("charge", d3.forceManyBody().strength(-1500))
        .force("center", d3.forceCenter(width / 2, height / 2));

    const link = zoomGroup.append("g")
        .attr("class", "links")
        .selectAll("line")
        .data(graph.links)
        .enter().append("line")
        .attr("stroke-width", 8)
        .attr("stroke", "#999");

    // Label edges with weights
    const linkLabels = zoomGroup.append("g")
        .attr("class", "link-labels")
        .selectAll("text")
        .data(graph.links)
        .enter().append("text")
        .attr("text-anchor", "middle")
        .attr("dy", -5)
        .attr("font-size", "12px")
        .attr("fill", "black")
        .text(d => `${d.weight.toFixed(3)} (${d.betw.toFixed(3)})`);

    const node = zoomGroup.append("g")
        .attr("class", "nodes")
        .selectAll("g")
        .data(graph.nodes)
        .enter().append("g");

    const circles = node.append("circle")
        .attr("r", 10)
        .attr("fill", function(d) {
            if (source_id_array.includes(d.id)) {
                return "green";
            } else if (sink_id_array.includes(d.id)) {
                return "red";
            } else {
                return "black";
            }
        })
        .on("click", function(event, d) {
            if (resistor_1 === d.id) {
                resistor_1 = null;
                d3.select(this).attr("fill", "black");
            } else if (resistor_2 === d.id) {
                resistor_2 = null;
                d3.select(this).attr("fill", "black");
            } else if (resistor_1 === null) {
                resistor_1 = d.id;
                d3.select(this).attr("fill", "yellow");
            } else if (resistor_2 === null) {
                resistor_2 = d.id;
                d3.select(this).attr("fill", "yellow");
            }

            updateCalculateButtonVisibility();

            console.log("Resistor 1 ID:", resistor_1, "Resistor 2 ID:", resistor_2);
        });

    const labels = node.append("text")
        .attr("x", 15)
        .attr("y", 5)
        .attr("fill", "black")
        .attr("font-size", "12px")
        .attr("text-anchor", "middle")
        .text(d => startingIndexValue === "1" ? d.id + 1 : d.id);

    simulation
        .nodes(graph.nodes)
        .on("tick", ticked);

    simulation.force("link")
        .links(graph.links);

    function ticked() {
        link
            .attr("x1", d => d.source.x)
            .attr("y1", d => d.source.y)
            .attr("x2", d => d.target.x)
            .attr("y2", d => d.target.y);

        linkLabels
            .attr("x", d => (d.source.x + d.target.x) / 2)
            .attr("y", d => (d.source.y + d.target.y) / 2);

        node
            .attr("transform", d => `translate(${d.x},${d.y})`);
    }

    function updateCalculateButtonVisibility() {
        const calculateButton = document.getElementById('calculate-button');
        if (resistor_1 !== null && resistor_2 !== null) {
            calculateButton.style.display = 'inline'; // Show button
        } else {
            calculateButton.style.display = 'none'; // Hide button
        }
    }

    // Handle Calculate button click
    document.getElementById('calculate-button').onclick = async function() {
        if (resistor_1 !== null && resistor_2 !== null) {
            const resist1Resid = resistor_1;
            const resist2Resid = resistor_2;

            console.log(resist1Resid, " is resist1");
            console.log(resist2Resid, " is resist2");

            const formData = new FormData();
            formData.append('resist1', resist1Resid);
            formData.append('resist2', resist2Resid);
            const response = axios.post('/api/calculate', formData, {
                headers: {
                    'Content-Type': 'application/json',
                },
            })
            
            const result = response.data;
            if (result.betweenness_score === 0) {
                alert("Please select 2 Resistors that have a direct connection");
            } else {
                document.getElementById('output').textContent = JSON.stringify(result, null, 2);
            }
        }
    };
    document.getElementById('refresh-button').onclick = async function() {
        resistor_1 = null;
        resistor_2 = null;
        drawGraph(graph);
        setupColorScaleAndEdges();
    };

}        

function displayRankedNodes(ranked_nodes, ranked_nodes2) {
    const container = d3.select('#ranked-nodes');
    container.selectAll("*").remove(); // Clear previous paths

    // Create a flex container to hold paths and paths2 side-by-side
    const flexContainer = container.append("div")
        .style("display", "flex"); // Use flexbox for horizontal alignment

    // Create a container for paths
    const pathsContainer = flexContainer.append("div")
        .html("<strong>Most Important Nodes Based off of Betweenness</strong>")
        .style("flex", "1");
        //.style("margin-right", "20px"); // Optional: Add space between paths and paths2


    ranked_nodes.forEach((node, index) => {
        // Create a div to hold the path information and button
        const pathDiv = pathsContainer.append("div")
            .style("display", "flex") // Use flexbox for horizontal alignment
            .style("align-items", "center") // Center items vertically
            .style("margin-bottom", "10px"); // Optional: Add space between rows
        
        pathDiv.append("div")
            .text(`${index + 1}: Node: ${startingIndexValue === "1" ? node.node + 1: node.node}, Frequency: ${node.frequency}`)
            .style("margin-right", "10px"); // Optional: Add space between text and button

        // Add a button to the div
        // pathDiv.append("button")
        //     .text("Highlight") // Change this text to whatever you want the button to display
        //     .on("click", () => {
        //         // Define what should happen when the button is clicked
        //         //alert(`Button for Path ${index + 1} clicked!`);
        //         //highlightPathEdges(path);
        //     });
    });

    const paths2Container = flexContainer.append("div")
        .html("<strong>Most Important Nodes Based off of Correlation</strong>")
        .style("flex", "1");

    ranked_nodes2.forEach((node, index) => {
        // Create a div to hold the path information and button
        const pathDiv = paths2Container.append("div")
            .style("display", "flex") // Use flexbox for horizontal alignment
            .style("align-items", "center") // Center items vertically
            .style("margin-bottom", "10px"); // Optional: Add space between rows
        
        pathDiv.append("div")
            .text(`${index + 1}: Node: ${startingIndexValue === "1" ? node.node + 1: node.node}, Frequency: ${node.frequency}`)
            .style("margin-right", "10px"); // Optional: Add space between text and button

        // Add a button to the div
        // pathDiv.append("button")
        //     .text("Highlight") // Change this text to whatever you want the button to display
        //     .on("click", () => {
        //         // Define what should happen when the button is clicked
        //         //alert(`Button for Path ${index + 1} clicked!`);
        //         //highlightPathEdges(path);
        //     });
    });
}

// Function to update the source node label
function updateSourceNodeLabel() {
    const sourceNodeLabel = document.getElementById('source-node-label');
    const sinkNodeLabel = document.getElementById('sink-node-label');
    const dirLabel = document.getElementById('directions-label');
    const refreshButton = document.getElementById('refresh-button');
    const downloadButton = document.getElementById('download-pdf');
    
    // Show the label
    sourceNodeLabel.style.display = 'inline';
    sinkNodeLabel.style.display = 'inline'; 
    dirLabel.style.display = 'inline';
    refreshButton.style.display = 'inline';
    downloadButton.style.display = 'inline';

}

function setupColorScaleAndEdges() {
    // Define the color scale
    const colorScale = d3.scaleLinear()
        .domain([0, large_bet]) // Using the maximum edge frequency
        .range(["lightgray", "red"]); // Low frequency -> black, High frequency -> red

    // Implement the colorEdges Function
    function colorEdges() {
        console.log("Color scale domain:", colorScale.domain());
        console.log("Color scale range:", colorScale.range());

        d3.selectAll(".links line") // Select all line elements in the links group
            .attr("stroke-width", 8)
            .attr("stroke", d => colorScale(d.betw)); // Set stroke color based on betweenness
    }
    colorEdges();
}

function drawColorScale() {
    const svg = d3.select("svg");
    
    // Remove any existing color scales
    svg.selectAll(".color-scale").remove();
    
    // Define dimensions and position for the color scale
    const width = 20;
    const height = 300;
    const x = svg.attr("width") - 60; // Positioning on the right side
    const y = 50;
    
    // Create a group for the color scale
    const colorScaleGroup = svg.append("g")
        .attr("class", "color-scale")
        .attr("transform", `translate(${x}, ${y})`);
    
    // Define the gradient
    const gradient = colorScaleGroup.append("defs")
        .append("linearGradient")
        .attr("id", "color-gradient")
        .attr("x1", "0%")
        .attr("y1", "100%")
        .attr("x2", "0%")
        .attr("y2", "0%");
    
    gradient.append("stop")
        .attr("offset", "0%")
        .attr("stop-color", "lightgray"); // Lightest
    
    gradient.append("stop")
        .attr("offset", "100%")
        .attr("stop-color", d3.interpolateReds(1)); // Darkest
    
    // Draw the rectangle filled with the gradient
    colorScaleGroup.append("rect")
        .attr("width", width)
        .attr("height", height)
        .style("fill", "url(#color-gradient)");
    
    // Define the scale for the color scale's axis
    const scale = d3.scaleLinear()
        .domain([0, large_bet])
        .range([height, 0]);
    
    // Define the axis for the color scale
    const axis = d3.axisRight(scale)
        .ticks(6); // Adjust the number of ticks as needed
    
    // Draw the axis
    colorScaleGroup.append("g")
        .attr("class", "axis")
        .attr("transform", `translate(${width}, 0)`)
        .call(axis);
}

function drawCorrelationMatrix(graph) {
    const margin = {top: 20, right: 0, bottom: 20, left: 50},
          width = 800 - margin.left - margin.right,
          height = 800 - margin.top - margin.bottom;

    const svg = d3.select("#correlation-svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);

    const nodes = graph.nodes.map(d => d.id);
    const nodeCount = nodes.length;
    const cellSize = width / nodeCount;

    // Define the color scale
    const colorScale = d3.scaleLinear()
        .domain([0, 1])  // Correlation values between 0 and 1
        .range(["red", "blue"]);

    // Create the matrix cells
    const cells = svg.selectAll("rect")
        .data(graph.links)
        .enter().append("rect")
        .attr("x", d => nodes.indexOf(d.source.id) * cellSize)
        .attr("y", d => nodes.indexOf(d.target.id) * cellSize)
        .attr("width", cellSize)
        .attr("height", cellSize)
        .attr("fill", d => colorScale(d.weight))
        .attr("stroke", "#ccc")
        .on("mouseover", function(event, d) {
            d3.select(this).attr("stroke", "black");
            const tooltip = d3.select("#correlation-svg").append("text")
                .attr("x", width + 10)
                .attr("y", height / 2)
                .attr("id", "tooltip")
                .attr("fill", "black")
                .text(`Nodes: ${d.source.id} & ${d.target.id}, Correlation: ${d.weight.toFixed(3)}`);
        })
        .on("mouseout", function() {
            d3.select(this).attr("stroke", "#ccc");
            d3.select("#tooltip").remove();
        });

    // Draw x and y axes (with limited labels)
    const xAxisScale = d3.scalePoint()
        .domain(nodes)
        .range([0, width]);

    const yAxisScale = d3.scalePoint()
        .domain(nodes)
        .range([0, height]);

    const xAxis = d3.axisBottom(xAxisScale)
        .tickValues(nodes.filter((d, i) => i % Math.ceil(nodeCount / 20) === 0));  // Show 20 evenly spaced labels

    const yAxis = d3.axisLeft(yAxisScale)
        .tickValues(nodes.filter((d, i) => i % Math.ceil(nodeCount / 20) === 0));  // Show 20 evenly spaced labels

    svg.append("g")
        .attr("transform", `translate(0,${height})`)
        .call(xAxis)
        .selectAll("text")
        .attr("transform", "rotate(-65)")
        .style("text-anchor", "end");

    svg.append("g")
        .call(yAxis);
}

document.getElementById('download-pdf').addEventListener('click', function() {
    // Select the contents of TopPaths and Histograms
    var topPathsContent = document.getElementById('TopPaths').innerHTML;
    var histogramsContent = document.getElementById('Histograms').innerHTML;

    // Create a container for the combined content
    var content = `
        <div>
            ${topPathsContent}
        </div>
        <div style="page-break-before: always;">
            ${histogramsContent}
        </div>
    `;

    // Create a hidden div element and add the content
    var element = document.createElement('div');
    element.innerHTML = content;
    document.body.appendChild(element);

    // Generate the PDF
    html2pdf().from(element).save('download.pdf').then(function() {
        // Remove the hidden element after generating the PDF
        document.body.removeChild(element);
    });
});

