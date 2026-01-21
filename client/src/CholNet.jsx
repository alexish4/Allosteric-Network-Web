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
  const EXAMPLE_PDB_URL = "7D93.pdb";

  const clrColorMap = useMemo(() => {
    if (!results || !Array.isArray(results)) return {};

    // stable order
    const sorted = [...results].sort(
      (a, b) =>
        (a.clr_residue_number ?? 0) - (b.clr_residue_number ?? 0) ||
        String(a.clr_chain_id ?? "").localeCompare(String(b.clr_chain_id ?? ""))
    );

    // HSL -> HEX helper (3Dmol likes hex)
    const hslToHex = (h, s, l) => {
      s /= 100;
      l /= 100;

      const k = (n) => (n + h / 30) % 12;
      const a = s * Math.min(l, 1 - l);
      const f = (n) =>
        l - a * Math.max(-1, Math.min(k(n) - 3, Math.min(9 - k(n), 1)));

      const toHex = (x) =>
        Math.round(255 * x)
          .toString(16)
          .padStart(2, "0");

      return `#${toHex(f(0))}${toHex(f(8))}${toHex(f(4))}`;
    };

    const map = {};
    const n = Math.max(sorted.length, 1);

    sorted.forEach((r, i) => {
      const chain = r?.clr_chain_id;
      const resi = r?.clr_residue_number;
      if (chain == null || resi == null) return;

      const hue = Math.round((i * 360) / n);
      map[`${chain}${resi}`] = hslToHex(hue, 80, 45); // hex color
    });

    return map;
  }, [results]);



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

  // Render whenever pdbText or results change (so CLR colors update after prediction)
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

    // Protein default
    viewer.setStyle({}, { cartoon: {} });

    // Default CLR style first (fallback)
    viewer.setStyle({ resn: "CLR" }, { stick: { color: "gold" } });

    // Color each CLR uniquely (hex colors)
    if (results && Array.isArray(results) && results.length > 0) {
      for (const r of results) {
        const chain = r?.clr_chain_id;
        const resi = r?.clr_residue_number;
        if (chain == null || resi == null) continue;

        const key = `${chain}${resi}`;
        const color = clrColorMap[key];
        if (!color) continue;

        viewer.setStyle(
          { resn: "CLR", chain: String(chain), resi: Number(resi) },
          { stick: { color: color } } // hex, e.g. "#aabbcc"
        );
      }
    }


    // Hover labels on CLR atoms
    viewer.setHoverable(
      { resn: "CLR" },
      true,
      function (atom) {
        if (atom.label) return;

        const chain = atom.chain || "";
        const resi = atom.resi != null ? atom.resi : "";
        const resn = atom.resn || "";
        const aname = atom.atom || atom.elem || "";

        atom.label = viewer.addLabel(`${resn} ${chain}${resi} • ${aname}`, {
          position: { x: atom.x, y: atom.y, z: atom.z },
          backgroundColor: "white",
          borderColor: "#333",
          borderThickness: 1,
          fontColor: "#111",
          fontSize: 12,
          inFront: true,
        });

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
  }, [pdbText, results, clrColorMap]); // <-- important: include results


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
        <div style={styles.description}>
          Please upload a PDB file containing both protein and docked cholesterol molecules.
          CholBindNet will provide a confidence ranking.
        </div>

      <div style={styles.uploadControls}>

        <div style={styles.exampleRow}>
          <a href={EXAMPLE_PDB_URL} download style={styles.exampleLink}>
            Download example PDB
          </a>
          <span style={styles.exampleHint}>
            (Use this to test the interface if you don’t have a file handy.)
          </span>
        </div>

        <form onSubmit={onSubmit} style={styles.formRow}>
          <label htmlFor="pdb_file" style={styles.label}>
            Select .pdb file:
          </label>
          <input
            type="file"
            id="pdb_file"
            accept=".pdb"
            required
            onChange={(e) => onFilePicked(e.target.files?.[0] || null)}
            style={styles.fileInput}
          />
          <button type="submit" style={styles.button} disabled={isLoading}>
            {isLoading ? "Analyzing..." : "Analyze"}
          </button>
        </form>

      </div>


        <div style={styles.viewerNote}>
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
          <h2>Results (CLR ligands)</h2>

          <div style={{ overflowX: "auto" }}>
            <table style={styles.table}>
              <thead>
                <tr>
                  <th style={styles.th}>CLR ID</th>
                  {modelOrder.map((m) => (
                    <th key={m} style={styles.th}>{m}</th>
                  ))}
                </tr>
              </thead>

              <tbody>
                {results.map((r, idx) => {
                  const clrKey = `${r.clr_chain_id}${r.clr_residue_number}`;
                  const clrColor = clrColorMap[clrKey] || "gold";

                  return (
                    <tr key={idx}>
                      <td style={styles.td}>
                        <span
                          title={`Color for CLR ${clrKey}`}
                          style={{
                            display: "inline-block",
                            width: 12,
                            height: 12,
                            borderRadius: 3,
                            background: clrColor,
                            border: "1px solid #999",
                            marginRight: 8,
                            verticalAlign: "middle",
                          }}
                        />
                        <strong>CLR {clrKey}</strong>
                      </td>

                      {modelOrder.map((modelName) => {
                        const data = r[modelName];
                        if (!data) {
                          return (
                            <td key={modelName} style={styles.tdMuted}>
                              —
                            </td>
                          );
                        }

                        const score = Number(data.mean_score);
                        const isHigh = score > 0.5;

                        return (
                          <td key={modelName} style={styles.td}>
                            <div style={{ fontWeight: 700, color: isHigh ? "green" : "#d9534f" }}>
                              {Number.isFinite(score) ? score.toFixed(4) : "—"}
                            </div>
                            <div style={{ fontSize: 12, color: "#555" }}>
                              {isHigh ? "Positive" : "Negative"}
                            </div>
                          </td>
                        );
                      })}
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>

          <div style={{ marginTop: 10, fontSize: 13, color: "#555" }}>
            Tip: each CLR row color matches the cholesterol color in the 3D viewer.
          </div>
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
    lineHeight: 1.7,
    fontSize: 18,
  },
  h1: { color: "#333", fontSize: 34, marginBottom: 10  },
  uploadSection: {
    background: "#f9f9f9",
    padding: 22,
    borderRadius: 10,
    marginBottom: 12,
    border: "1px solid #ddd",
  },
  description: {
    fontSize: 18,
    color: "#333",
    marginBottom: 10,
  },
  uploadControls: {
    display: "flex",
    flexDirection: "column",
    alignItems: "center",     // ⭐ centers the whole block
    gap: 14,
    marginTop: 10,
  },

  exampleRow: {
    display: "flex",
    alignItems: "center",
    justifyContent: "center", // ⭐ centers example row
    gap: 10,
    flexWrap: "wrap",
    textAlign: "center",
  },

  exampleLink: {
    display: "inline-block",
    padding: "8px 12px",
    borderRadius: 8,
    border: "1px solid #ccc",
    background: "white",
    color: "#111",
    textDecoration: "none",
    fontWeight: 700,
  },

  exampleHint: {
    color: "#555",
    fontSize: 14,
  },

  formRow: {
    display: "flex",
    alignItems: "center",
    justifyContent: "center", // ⭐ centers file chooser + button
    gap: 12,
    flexWrap: "wrap",
  },

  label: {
    fontWeight: 700,
    fontSize: 18,
  },

  fileInput: {
    fontSize: 16,
  },
  error: { color: "red", background: "#ffe6e6", padding: 10, borderRadius: 4, marginBottom: 12 },
  button: {
    padding: "10px 18px",
    cursor: "pointer",
    fontSize: 16,
    fontWeight: 700,
    borderRadius: 8,
    border: "1px solid #bbb",
    background: "white",
  },
  viewerNote: {
    marginTop: 12,
    fontSize: 16,
    color: "#444",
  },
  resultsSection: { marginTop: 20 },
  fileBlock: { marginTop: 20, borderBottom: "2px solid #eee", paddingBottom: 20 },
  modelCard: { border: "1px solid #ccc", padding: 15, marginBottom: 15, borderRadius: 6 },
  modelName: { fontWeight: "bold", fontSize: "1.1em" },
  score: { fontSize: "1.3em", fontWeight: "bold" },
  table: {
    width: "100%",
    borderCollapse: "collapse",
    marginTop: 10,
    background: "white",
    border: "1px solid #ddd",
    borderRadius: 8,
  },
  th: {
    textAlign: "left",
    padding: "10px 12px",
    borderBottom: "1px solid #ddd",
    background: "#f7f7f7",
    fontWeight: 700,
    whiteSpace: "nowrap",
  },
  td: {
    padding: "10px 12px",
    borderBottom: "1px solid #eee",
    verticalAlign: "top",
    whiteSpace: "nowrap",
  },
  tdMuted: {
    padding: "10px 12px",
    borderBottom: "1px solid #eee",
    color: "#999",
    whiteSpace: "nowrap",
  },

};
