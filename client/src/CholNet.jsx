import React, { useEffect, useMemo, useRef, useState } from "react";
import axios from "axios";

export default function CholNet() {
  const [file, setFile] = useState(null);
  const [pdbText, setPdbText] = useState("");
  const [error, setError] = useState("");
  const [results, setResults] = useState(null);
  const [isLoading, setIsLoading] = useState(false);

  const viewerDivRef = useRef(null);
  const viewerRef = useRef(null);
  const $3DmolRef = useRef(null);

  const modelOrder = useMemo(() => ["GAT", "GCN", "GNN"], []);

  // Load 3Dmol once
  useEffect(() => {
    let mounted = true;
    import("3dmol/build/3Dmol.js")
      .then(($3Dmol) => {
        if (!mounted) return;
        $3DmolRef.current = $3Dmol;
      })
      .catch((e) => {
        console.error(e);
        setError("Failed to load 3D viewer (3Dmol).");
      });
    return () => {
      mounted = false;
    };
  }, []);

  // Read the PDB locally on file select (so it renders instantly)
  const onFilePicked = (f) => {
    setError("");
    setResults(null);
    setFile(f || null);

    if (!f) return;
    if (!f.name.endsWith(".pdb")) {
      setError("Invalid file type. Please upload a .pdb file.");
      setPdbText("");
      return;
    }

    const reader = new FileReader();
    reader.onload = () => setPdbText(String(reader.result || ""));
    reader.onerror = () => setError("Failed to read the .pdb file.");
    reader.readAsText(f);
  };

  // Render whenever pdbText changes
  useEffect(() => {
    const $3Dmol = $3DmolRef.current;
    const el = viewerDivRef.current;
    if (!$3Dmol || !el || !pdbText) return;

    // Create viewer once, then reuse it
    if (!viewerRef.current) {
      viewerRef.current = $3Dmol.createViewer(el, { backgroundColor: "white" });
    }
    const viewer = viewerRef.current;

    viewer.clear();
    viewer.addModel(pdbText, "pdb");

    // Protein / everything default: cartoon (you can tweak later)
    viewer.setStyle({}, { cartoon: {} });

    // Cholesterol as bond/stick representation: resn "CLR"
    // (This is the key line you asked for.)
    viewer.setStyle(
      { resn: "CLR" },
      { stick: { color: "gold" } }   // or "orange", "cyan", "#FFD700"
    );

    viewer.setHoverable(
      { resn: "CLR" },          // only cholesterol atoms get hover
      true,
      function (atom) {
        if (atom.label) return;

        const chain = atom.chain || "";
        const resi = atom.resi != null ? atom.resi : "";
        const resn = atom.resn || "";
        const aname = atom.atom || atom.elem || "";

        atom.label = viewer.addLabel(
          `${resn} ${chain}${resi} • ${aname}`,
          {
            position: { x: atom.x, y: atom.y, z: atom.z },
            backgroundColor: "white",
            borderColor: "#333",
            borderThickness: 1,
            fontColor: "#111",
            fontSize: 12,
            inFront: true,
          }
        );

        viewer.render();
      },
      function (atom) {
        if (atom.label) {
          viewer.removeLabel(atom.label);
          delete atom.label;
          viewer.render();
        }
      }
    );

    viewer.zoomTo();
    viewer.render();
    viewer.resize();
  }, [pdbText]);

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
      <h1 style={styles.h1}>CholNet Interface</h1>

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
            onChange={(e) => onFilePicked(e.target.files?.[0] || null)}
          />{" "}
          <button type="submit" style={styles.button} disabled={isLoading}>
            {isLoading ? "Analyzing..." : "Analyze"}
          </button>
        </form>

        <div style={{ marginTop: 12, fontSize: 14, color: "#444" }}>
          Viewer: protein shown as cartoon, cholesterol (<code>CLR</code>) shown as bonds/sticks.
        </div>
      </div>

      {/* 3Dmol viewport (full-bleed / breaks out of centered containers) */}
      <div
        style={{
          width: "100vw",
          marginLeft: "calc(50% - 50vw)",
          marginRight: "calc(50% - 50vw)",
          padding: "0 20px", // optional breathing room from screen edge
        }}
      >
        <div
          ref={viewerDivRef}
          style={{
            width: "100%",
            height: 520,
            border: "1px solid #ddd",
            borderRadius: 8,
            background: "white",
            position: "relative",
            overflow: "hidden",
          }}
        />
      </div>


      {results && (
        <div style={styles.resultsSection}>
          <h2>Results (per CLR ligand)</h2>

          {results.map((r, idx) => (
            <div key={idx} style={styles.fileBlock}>
              <h3>
                File: {r.filename} — CLR {r.clr_chain_id}{r.clr_residue_number}
              </h3>

              {modelOrder.map((modelName) => {
                const data = r[modelName];
                if (!data) {
                  return (
                    <div key={modelName} style={styles.modelCard}>
                      <div style={styles.modelName}>{modelName} Model</div>
                      <div style={{ color: "gray" }}>
                        Result unavailable (processing failed or skipped)
                      </div>
                    </div>
                  );
                }

                const score = data.mean_score;
                const isHigh = score > 0.5;

                return (
                  <div key={modelName} style={styles.modelCard}>
                    <div style={styles.modelName}>{modelName} Model</div>
                    <div style={{ ...styles.score, color: isHigh ? "green" : "#d9534f" }}>
                      Probability: {Number(score).toFixed(4)}
                    </div>
                    <div>
                      Prediction:{" "}
                      {isHigh ? <strong>Positive (Binding)</strong> : <span>Negative</span>}
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
    maxWidth: "95vw",   // nearly full screen
    margin: "0 auto",
    padding: 20,
    lineHeight: 1.6,
  },
  h1: { color: "#333" },
  uploadSection: {
    background: "#f9f9f9",
    padding: 20,
    borderRadius: 8,
    marginBottom: 12,
    border: "1px solid #ddd",
  },
  error: { color: "red", background: "#ffe6e6", padding: 10, borderRadius: 4, marginBottom: 12 },
  button: { padding: "5px 15px", cursor: "pointer" },
  resultsSection: { marginTop: 20 },
  fileBlock: { marginTop: 20, borderBottom: "2px solid #eee", paddingBottom: 20 },
  modelCard: { border: "1px solid #ccc", padding: 15, marginBottom: 15, borderRadius: 6 },
  modelName: { fontWeight: "bold", fontSize: "1.1em" },
  score: { fontSize: "1.3em", fontWeight: "bold" },
};
