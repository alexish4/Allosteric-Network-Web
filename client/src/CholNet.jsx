import React, { useMemo, useState } from "react";
import axios from "axios";

export default function CholNet() {
  const [file, setFile] = useState(null);
  const [error, setError] = useState("");
  const [results, setResults] = useState(null);
  const [isLoading, setIsLoading] = useState(false);

  const modelOrder = useMemo(() => ["GAT", "GCN", "GNN"], []);

  const onSubmit = async (e) => {
    e.preventDefault();
    setError("");
    setResults(null);

    if (!file) {
      setError("Please choose a .pdb file.");
      return;
    }
    if (!file.name.endsWith(".pdb")) {
      setError("Invalid file type. Please upload a .pdb file.");
      return;
    }

    const formData = new FormData();
    formData.append("pdb_file", file);

    setIsLoading(true);
    try {
      const res = await axios.post("/api/cholnet", formData, {
        headers: { "Content-Type": "multipart/form-data" },
      });

      if (res.data?.status === "success") {
        setResults(res.data.results);
      } else {
        setError(res.data?.message || "Processing error.");
      }
    } catch (err) {
      const msg =
        err?.response?.data?.message ||
        err?.message ||
        "An error occurred while processing the file.";
      setError(msg);
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <div style={styles.body}>
      <h1 style={styles.h1}>CholNet Inference</h1>
      <p>
        Upload a PDB file to predict cholesterol binding sites using GAT, GCN, and GNN models.
        <br /><br />
        <strong>File requirements:</strong> The PDB must contain cholesterol molecules annotated as
        <code> CLR </code> residues using standard <code>HETATM</code> records.
        <br />
        Example:
        <br />
        <code style={{ display: "block", marginTop: "6px" }}>
            HETATM19758  C1  CLR A1001     113.576 122.575 109.835  1.00 67.75           C
        </code>
      </p>

      {error && (
        <div style={styles.error}>
          <strong>Error:</strong> {error}
        </div>
      )}

      <div style={styles.uploadSection}>
        <form onSubmit={onSubmit}>
          <label htmlFor="pdb_file">Select .pdb file:</label>{" "}
          <input
            type="file"
            id="pdb_file"
            accept=".pdb"
            required
            onChange={(e) => setFile(e.target.files?.[0] || null)}
          />{" "}
          <button type="submit" style={styles.button} disabled={isLoading}>
            {isLoading ? "Analyzing..." : "Analyze"}
          </button>
        </form>
      </div>

      {results && (
        <div style={styles.resultsSection}>
          <h2>Results</h2>

          {results.map((result, idx) => (
            <div key={idx} style={styles.fileBlock}>
              <h3>File: {result.filename}</h3>

              {modelOrder.map((modelName) => {
                const data = result[modelName];
                if (!data) {
                  return (
                    <div key={modelName} style={styles.modelCard}>
                      <div style={styles.modelName}>{modelName} Model</div>
                      <div style={{ color: "gray" }}>
                        Result unavailable (processing failed or file skipped)
                      </div>
                    </div>
                  );
                }

                const score = data.mean_score;
                const isHigh = score > 0.5;

                // backend might send experiment_breakdown as object/map
                const breakdownVals = data.experiment_breakdown
                  ? Object.values(data.experiment_breakdown).map((v) =>
                      typeof v === "number" ? v.toFixed(4) : String(v)
                    )
                  : [];

                return (
                  <div key={modelName} style={styles.modelCard}>
                    <div style={styles.modelName}>{modelName} Model</div>

                    <div style={{ ...styles.score, color: isHigh ? "green" : "#d9534f" }}>
                      Probability: {Number(score).toFixed(4)}
                    </div>

                    <div>
                      Prediction:{" "}
                      {isHigh ? (
                        <strong>Positive (Binding)</strong>
                      ) : (
                        <span>Negative (Non-Binding)</span>
                      )}
                    </div>

                    <div style={styles.details}>
                      Based on average of ensemble models.
                      {breakdownVals.length > 0 && (
                        <>
                          <br />
                          Experiment scores: {breakdownVals.join(", ")}
                        </>
                      )}
                    </div>
                  </div>
                );
              })}
            </div>
          ))}
        </div>
      )}
    </div>
  );
}

const styles = {
  body: {
    fontFamily: "sans-serif",
    maxWidth: 800,
    margin: "0 auto",
    padding: 20,
    lineHeight: 1.6,
  },
  h1: { color: "#333" },
  uploadSection: {
    background: "#f9f9f9",
    padding: 20,
    borderRadius: 8,
    marginBottom: 20,
    border: "1px solid #ddd",
  },
  error: {
    color: "red",
    background: "#ffe6e6",
    padding: 10,
    borderRadius: 4,
    marginBottom: 20,
  },
  button: { padding: "5px 15px", cursor: "pointer" },
  resultsSection: { marginTop: 30 },
  fileBlock: {
    marginTop: 20,
    borderBottom: "2px solid #eee",
    paddingBottom: 20,
  },
  modelCard: {
    border: "1px solid #ccc",
    padding: 15,
    marginBottom: 15,
    borderRadius: 6,
  },
  modelName: { fontWeight: "bold", fontSize: "1.2em", color: "#0056b3" },
  score: { fontSize: "1.5em", fontWeight: "bold" },
  details: { fontSize: "0.9em", color: "#666", marginTop: 5 },
};
