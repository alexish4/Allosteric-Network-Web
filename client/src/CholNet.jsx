import React, { useEffect, useMemo, useRef, useState } from "react";
import axios from "axios";

export default function CholNet() {
  const [file, setFile] = useState(null);
  const [pdbText, setPdbText] = useState("");
  const [error, setError] = useState("");
  const [results, setResults] = useState(null);
  const [isLoading, setIsLoading] = useState(false);

  const [files, setFiles] = useState([]); // batch PDBs
  const [batchResults, setBatchResults] = useState(null);
  const [batchInterpretation, setBatchInterpretation] = useState(null);
  const [selectedClr, setSelectedClr] = useState(null);
  const [selectedModel, setSelectedModel] = useState("GNN");
  const [interpretationSchemaVersion, setInterpretationSchemaVersion] =
    useState(null);
  const isBatchMode = files.length > 0;

  const viewerDivRef = useRef(null);
  const viewerRef = useRef(null);
  const $3DmolRef = useRef(null);

  const modelOrder = useMemo(() => ["GNN", "GAT", "GCN"], []);

  // Average spy percentiles across Experiments 1–5 for each model.
  const SCORE_THRESHOLDS = useMemo(
    () => ({
      GNN: {
        p25: 0.6876584,
        p50: 0.7778270,
        p75: 0.8349516,
      },
      GAT: {
        p25: 0.6400374,
        p50: 0.8283712,
        p75: 0.9087056,
      },
      GCN: {
        p25: 0.6724110,
        p50: 0.8102502,
        p75: 0.8893382,
      },
    }),
    []
  );

  const getScoreLabel = (modelName, value) => {
    const score = Number(value);
    const thresholds = SCORE_THRESHOLDS[modelName];

    if (!Number.isFinite(score) || !thresholds) return "Unavailable";
    if (score < thresholds.p25) return "Negative";
    if (score < thresholds.p50) return "PseudoNegative";
    if (score < thresholds.p75) return "PseudoPositive";
    return "Positive";
  };

  const getLabelColor = (label) => {
    switch (label) {
      case "Negative":
        return "#d9534f";
      case "PseudoNegative":
        return "#d97706";
      case "PseudoPositive":
        return "#2563eb";
      case "Positive":
        return "#15803d";
      default:
        return "#777";
    }
  };
  
  const EXAMPLE_PDB_URL = "7D93.pdb";
  const ABSTRACT_FIGURE_URL = "CholBindAbstractFigure.png";
  const PUBLICATION_URL = "https://doi.org/10.1038/s42004-026-02064-w";

  const getGnnScore = (r) => {
    const v = Number(r?.GNN?.mean_score);
    return Number.isFinite(v) ? v : -Infinity;
  };

  const sortedResults = useMemo(() => {
    if (!Array.isArray(results)) return [];
    return [...results].sort((a, b) => getGnnScore(b) - getGnnScore(a));
  }, [results]);

  const rankedResults = useMemo(() => {
    return sortedResults.map((r, idx) => ({
      ...r,
      rank: idx + 1,
    }));
  }, [sortedResults]);

  const selectedResult = useMemo(() => {
    if (!selectedClr) return null;
    return (
      rankedResults.find(
        (result) =>
          `${result?.clr_chain_id ?? ""}${result?.clr_residue_number ?? ""}` ===
          selectedClr.key
      ) || null
    );
  }, [rankedResults, selectedClr]);

  const selectedInterpretation =
    selectedResult?.[selectedModel]?.interpretation || null;
  const selectedModelResult = selectedResult?.[selectedModel] || null;

  const clrColorMap = useMemo(() => {
    if (!Array.isArray(results)) return {};

    const sorted = [...results].sort(
      (a, b) =>
        (a.clr_residue_number ?? 0) - (b.clr_residue_number ?? 0) ||
        String(a.clr_chain_id ?? "").localeCompare(String(b.clr_chain_id ?? ""))
    );

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
      map[`${chain}${resi}`] = hslToHex(hue, 80, 45);
    });

    return map;
  }, [results]);

  const chooseClr = (r) => {
    const chain = String(r?.clr_chain_id ?? "");
    const resi = Number(r?.clr_residue_number);

    if (!Number.isFinite(resi)) return;

    const key = `${chain}${resi}`;

    // Clicking the selected CLR again clears the selection.
    if (selectedClr?.key === key) {
      setSelectedClr(null);
      return;
    }

    const firstAvailableModel =
      modelOrder.find((modelName) => r?.[modelName]?.interpretation) ||
      modelOrder.find((modelName) => r?.[modelName]);
    setSelectedModel(firstAvailableModel || "GNN");
    setSelectedClr({ chain, resi, key });
  };

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

  // Single-file picker: also clears batch mode
  const onFilePicked = (f) => {
    setError("");
    setResults(null);
    setBatchResults(null);
    setBatchInterpretation(null);
    setSelectedClr(null);
    setInterpretationSchemaVersion(null);
    setFiles([]); // leave batch mode

    setFile(f || null);
    if (!f) {
      setPdbText("");
      return;
    }

    if (!f.name.toLowerCase().endsWith(".pdb")) {
      setError("Invalid file type. Please upload a .pdb file.");
      setPdbText("");
      return;
    }

    const reader = new FileReader();
    reader.onload = () => setPdbText(String(reader.result || ""));
    reader.onerror = () => setError("Failed to read the .pdb file.");
    reader.readAsText(f);
  };

  // Viewer render (single mode only, because we clear pdbText in batch mode and we hide viewer)
  useEffect(() => {
    const $3Dmol = $3DmolRef.current;
    const el = viewerDivRef.current;
    if (!$3Dmol || !el || !pdbText) return;

    if (!viewerRef.current) {
      viewerRef.current = $3Dmol.createViewer(el, { backgroundColor: "white" });
    }
    const viewer = viewerRef.current;

    viewer.clear();
    viewer.addModel(pdbText, "pdb");

    viewer.setStyle({}, { cartoon: {} });
    viewer.setStyle(
      { resn: "CLR" },
      { stick: { color: "gold", radius: 0.2 } }
    );

    if (Array.isArray(results) && results.length > 0) {
      for (const r of results) {
        const chain = r?.clr_chain_id;
        const resi = r?.clr_residue_number;
        if (chain == null || resi == null) continue;

        const key = `${chain}${resi}`;
        const color = clrColorMap[key];
        if (!color) continue;

        viewer.setStyle(
          { resn: "CLR", chain: String(chain), resi: Number(resi) },
          { stick: { color, radius: 0.22 } }
        );
      }
    }

    // When a CLR is selected from the table, dim every CLR and strongly
    // highlight only the selected chain/residue combination.
    if (selectedClr) {
      viewer.setStyle(
        { resn: "CLR" },
        { stick: { color: "#c7c7c7", radius: 0.14 } }
      );

      const selectedClrSel = {
        resn: "CLR",
        chain: selectedClr.chain,
        resi: selectedClr.resi,
      };

      viewer.setStyle(selectedClrSel, {
        stick: { color: "#ff00ff", radius: 0.38 },
        sphere: { color: "#ff00ff", scale: 0.28 },
      });

      // Interpretation display: orange sticks mark important residues and red
      // spheres mark the exact important atoms for the selected architecture.
      const importantResidues = Array.isArray(
        selectedInterpretation?.top_residues
      )
        ? selectedInterpretation.top_residues
        : [];
      const importantAtoms = Array.isArray(selectedInterpretation?.top_atoms)
        ? selectedInterpretation.top_atoms
        : [];

      for (const residue of importantResidues) {
        const selection = {
          chain: String(residue?.chain_id ?? ""),
          resi: Number(residue?.residue_number),
        };
        if (!Number.isFinite(selection.resi)) continue;

        if (typeof viewer.addStyle === "function") {
          viewer.addStyle(selection, {
            stick: { color: "#f59e0b", radius: 0.18 },
          });
        }
      }

      for (const atom of importantAtoms) {
        const serial = Number(atom?.atom_serial);
        if (!Number.isFinite(serial)) continue;
        const selection = { serial };

        if (typeof viewer.addStyle === "function") {
          viewer.addStyle(selection, {
            stick: { color: "#dc2626", radius: 0.23 },
            sphere: { color: "#dc2626", scale: 0.34 },
          });
        }
      }
    }

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

    const cholSel = { resn: ["CLR", "CHL"] };
    const nearProteinSel = {
      and: [{ protein: true }, { within: { distance: 5.0, sel: cholSel } }],
    };

    if (typeof viewer.addStyle === "function") {
      viewer.addStyle(nearProteinSel, { sphere: { opacity: 0.25 } });
    } else {
      viewer.setStyle(nearProteinSel, {
        cartoon: { color: "spectrum" },
        vdw: { opacity: 0.25 },
      });
    }

    if (selectedClr) {
      viewer.zoomTo({
        resn: "CLR",
        chain: selectedClr.chain,
        resi: selectedClr.resi,
      });
    } else {
      viewer.zoomTo();
    }

    viewer.render();
    viewer.resize();
  }, [
    pdbText,
    results,
    clrColorMap,
    selectedClr,
    selectedModel,
    selectedInterpretation,
  ]);

  const downloadBatchCsv = () => {
    if (!Array.isArray(batchResults)) return;

    const rows = [];

    const summaryText = (modelName, categoryKey, labelKey) => {
      const records =
        batchInterpretation?.models?.[modelName]?.[categoryKey] || [];
      return records
        .map(
          (record) =>
            `${record[labelKey]} (${Number(record.percent || 0).toFixed(1)}%)`
        )
        .join("; ");
    };

    for (const item of batchResults) {
      const fname = item.filename || "";
      const resultsArr = Array.isArray(item.results) ? item.results : [];

      const ranked = [...resultsArr]
        .sort((a, b) => getGnnScore(b) - getGnnScore(a))
        .map((r, idx) => ({
          ...r,
          rank: idx + 1,
        }));

      for (const r of ranked) {
        const chain = r?.clr_chain_id ?? "";
        const resi = r?.clr_residue_number ?? "";
        const clr_id = `CLR ${chain}${resi}`;

        const gnnScore = r?.GNN?.mean_score;
        const gatScore = r?.GAT?.mean_score;
        const gcnScore = r?.GCN?.mean_score;

        rows.push({
          filename: fname,
          rank: r.rank,
          clr_id,
          GNN: gnnScore ?? "",
          GNN_label: getScoreLabel("GNN", gnnScore),
          GAT: gatScore ?? "",
          GAT_label: getScoreLabel("GAT", gatScore),
          GCN: gcnScore ?? "",
          GCN_label: getScoreLabel("GCN", gcnScore),
          GNN_top_atom_types: summaryText("GNN", "top_atom_types", "atom_type"),
          GNN_top_residue_types: summaryText(
            "GNN",
            "top_residue_types",
            "residue_name"
          ),
          GAT_top_atom_types: summaryText("GAT", "top_atom_types", "atom_type"),
          GAT_top_residue_types: summaryText(
            "GAT",
            "top_residue_types",
            "residue_name"
          ),
          GCN_top_atom_types: summaryText("GCN", "top_atom_types", "atom_type"),
          GCN_top_residue_types: summaryText(
            "GCN",
            "top_residue_types",
            "residue_name"
          ),
        });
      }
    }

    const headers = [
      "filename",
      "rank",
      "clr_id",
      "GNN",
      "GNN_label",
      "GAT",
      "GAT_label",
      "GCN",
      "GCN_label",
      "GNN_top_atom_types",
      "GNN_top_residue_types",
      "GAT_top_atom_types",
      "GAT_top_residue_types",
      "GCN_top_atom_types",
      "GCN_top_residue_types",
    ];
    const esc = (v) => {
      const s = String(v ?? "");
      return /[",\n]/.test(s) ? `"${s.replaceAll('"', '""')}"` : s;
    };

    const csv = [
      headers.join(","),
      ...rows.map((row) => headers.map((h) => esc(row[h])).join(",")),
    ].join("\n");

    const blob = new Blob([csv], { type: "text/csv;charset=utf-8;" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "cholnet_batch_results.csv";
    a.click();
    URL.revokeObjectURL(url);
  };

  const runSingle = async () => {
    setError("");
    setResults(null);
    setBatchResults(null);
    setBatchInterpretation(null);
    setSelectedClr(null);
    setInterpretationSchemaVersion(null);

    if (!file || !file.name.toLowerCase().endsWith(".pdb")) {
      setError("Please choose a .pdb file.");
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
        const returnedResults = Array.isArray(res.data.results)
          ? res.data.results
          : [];
        setResults(returnedResults);
        setInterpretationSchemaVersion(
          res.data.interpretation_schema_version ?? null
        );

        const highestRanked = [...returnedResults].sort(
          (a, b) => getGnnScore(b) - getGnnScore(a)
        )[0];
        if (highestRanked) {
          const chain = String(highestRanked.clr_chain_id ?? "");
          const resi = Number(highestRanked.clr_residue_number);
          const firstAvailableModel =
            modelOrder.find(
              (modelName) => highestRanked?.[modelName]?.interpretation
            ) || modelOrder.find((modelName) => highestRanked?.[modelName]);
          if (Number.isFinite(resi)) {
            setSelectedClr({ chain, resi, key: `${chain}${resi}` });
          }
          setSelectedModel(firstAvailableModel || "GNN");
        }
      } else setError(res.data?.message || "Processing error.");
    } catch (err) {
      setError(err?.response?.data?.message || err?.message || "Request failed.");
    } finally {
      setIsLoading(false);
    }
  };

  const runBatch = async () => {
    setError("");
    setResults(null);
    setBatchResults(null);
    setBatchInterpretation(null);
    setSelectedClr(null);
    setInterpretationSchemaVersion(null);

    if (!files || files.length === 0) {
      setError("Please select a folder (batch) with .pdb files.");
      return;
    }

    const formData = new FormData();
    for (const f of files) formData.append("pdb_files", f);

    setIsLoading(true);
    try {
      const res = await axios.post("/api/cholnet/batch", formData, {
        headers: { "Content-Type": "multipart/form-data" },
      });

      if (res.data?.status === "success") {
        setBatchResults(res.data.items);
        setBatchInterpretation(res.data.interpretation_summary || null);
      } else setError(res.data?.message || "Batch processing error.");
    } catch (err) {
      setError(err?.response?.data?.message || err?.message || "Batch request failed.");
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <div style={styles.body}>
      <h1 style={styles.h1}>CholBindNet Interface</h1>

      <div style={styles.publicationSection}>
        <a
          href={PUBLICATION_URL}
          target="_blank"
          rel="noopener noreferrer"
          style={styles.abstractFigureLink}
          title="View the CholBindNet publication"
        >
          <img
            src={ABSTRACT_FIGURE_URL}
            alt="CholBindNet graphical abstract"
            style={styles.abstractFigure}
          />
        </a>

        <div style={styles.publicationText}>
          <strong>CholBindNet</strong> is an interpretable neural-network framework
          for classifying cholesterol-binding sites in transmembrane proteins.
        </div>

        <a
          href={PUBLICATION_URL}
          target="_blank"
          rel="noopener noreferrer"
          style={styles.publicationLink}
        >
          Read the CholBindNet publication
        </a>
      </div>

      {error && (
        <div style={styles.error}>
          <strong>Error:</strong> {error}
        </div>
      )}

      <div style={styles.uploadSection}>
        <div style={styles.description}>
          Upload a PDB file containing both protein and docked cholesterol molecules.
          CholBindNet will provide a confidence ranking.
        </div>

        <div style={styles.uploadControls}>
          <div style={styles.exampleRow}>
            <a href={EXAMPLE_PDB_URL} download style={styles.exampleLink}>
              Download example PDB
            </a>
            <span style={styles.exampleHint}>(Use this to test the interface.)</span>
          </div>

          {/* ---------- Single mode ---------- */}
          <div style={styles.formRow}>
            <label htmlFor="pdb_file" style={styles.label}>
              Single PDB:
            </label>
            <input
              type="file"
              id="pdb_file"
              accept=".pdb"
              onChange={(e) => onFilePicked(e.target.files?.[0] || null)}
              style={styles.fileInput}
            />
            <button
              type="button"
              style={styles.button}
              onClick={runSingle}
              disabled={isLoading || isBatchMode}
              title={isBatchMode ? "Clear batch selection to run single-file" : ""}
            >
              {isLoading && !isBatchMode ? "Analyzing..." : "Analyze Single"}
            </button>
          </div>

          {/* ---------- Batch mode ---------- */}
          <div style={styles.formRow}>
            <label htmlFor="pdb_folder" style={styles.label}>
              Folder (batch):
            </label>

            <input
              type="file"
              id="pdb_folder"
              accept=".pdb"
              multiple
              webkitdirectory="true"
              directory="true"
              onChange={(e) => {
                setError("");
                setResults(null);
                setBatchResults(null);
                setBatchInterpretation(null);
                setSelectedClr(null);
                setInterpretationSchemaVersion(null);

                const picked = Array.from(e.target.files || []).filter((f) =>
                  f.name.toLowerCase().endsWith(".pdb")
                );

                setFiles(picked);
                setFile(null);
                setPdbText(""); // ensures viewer stays single-mode only
              }}
              style={styles.fileInput}
            />

            <button
              type="button"
              style={styles.button}
              onClick={runBatch}
              disabled={isLoading || files.length === 0}
            >
              {isLoading && isBatchMode ? "Running..." : "Run Batch"}
            </button>

            {files.length > 0 && (
              <div style={{ fontSize: 13, color: "#555" }}>
                {files.length} PDBs selected
                <button
                  type="button"
                  style={{ ...styles.button, marginLeft: 10, padding: "6px 10px" }}
                  onClick={() => {
                    setFiles([]);
                    setBatchResults(null);
                    setBatchInterpretation(null);
                    setSelectedClr(null);
                    setInterpretationSchemaVersion(null);
                  }}
                >
                  Clear
                </button>
              </div>
            )}
          </div>

          {!isBatchMode && (
            <div style={styles.viewerNote}>
              Viewer: protein shown as cartoon, cholesterol (<code>CLR</code>) shown as bonds/sticks.
            </div>
          )}
        </div>
      </div>

      {/* ---------- 3D viewer: SINGLE ONLY ---------- */}
      {!isBatchMode && (
        <div style={{ marginTop: 12 }}>
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
      )}

      {/* ---------- Single results ---------- */}
      {Array.isArray(results) && (
        <div style={styles.resultsSection}>
          <div style={styles.resultsHeader}>
            <h2 style={{ margin: 0 }}>Results (CLR ligands)</h2>

            {selectedClr && (
              <div style={styles.selectedClrControls}>
                <span>
                  Highlighting <strong>CLR {selectedClr.key}</strong>
                </span>
                <button
                  type="button"
                  style={{ ...styles.button, padding: "6px 10px" }}
                  onClick={() => setSelectedClr(null)}
                >
                  Show all CLRs
                </button>
              </div>
            )}
          </div>

          <div style={styles.thresholdLegend}>
            <div style={styles.thresholdTitle}>
              Model-specific thresholds averaged across five spy experiments
            </div>

            <div style={{ overflowX: "auto" }}>
              <table style={styles.thresholdTable}>
                <thead>
                  <tr>
                    <th style={styles.thresholdTh}>Model</th>
                    <th style={styles.thresholdTh}>Negative</th>
                    <th style={styles.thresholdTh}>PseudoNegative</th>
                    <th style={styles.thresholdTh}>PseudoPositive</th>
                    <th style={styles.thresholdTh}>Positive</th>
                  </tr>
                </thead>
                <tbody>
                  {modelOrder.map((modelName) => {
                    const t = SCORE_THRESHOLDS[modelName];

                    return (
                      <tr key={modelName}>
                        <td style={styles.thresholdTd}>
                          <strong>{modelName}</strong>
                        </td>
                        <td style={styles.thresholdTd}>
                          score &lt; {t.p25.toFixed(6)}
                        </td>
                        <td style={styles.thresholdTd}>
                          {t.p25.toFixed(6)} ≤ score &lt; {t.p50.toFixed(6)}
                        </td>
                        <td style={styles.thresholdTd}>
                          {t.p50.toFixed(6)} ≤ score &lt; {t.p75.toFixed(6)}
                        </td>
                        <td style={styles.thresholdTd}>
                          score ≥ {t.p75.toFixed(6)}
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>
          </div>

          <div style={{ overflowX: "auto" }}>
            <table style={styles.table}>
              <thead>
                <tr>
                  <th style={styles.th}>Rank</th>
                  <th style={styles.th}>CLR ID</th>
                  {modelOrder.map((m) => (
                    <th key={m} style={styles.th}>
                      {m}
                    </th>
                  ))}
                </tr>
              </thead>

              <tbody>
                {rankedResults.map((r, idx) => {
                  const clrKey = `${r.clr_chain_id}${r.clr_residue_number}`;
                  const clrColor = clrColorMap[clrKey] || "gold";

                  const isSelected = selectedClr?.key === clrKey;

                  return (
                    <tr
                      key={idx}
                      onClick={() => chooseClr(r)}
                      onKeyDown={(e) => {
                        if (e.key === "Enter" || e.key === " ") {
                          e.preventDefault();
                          chooseClr(r);
                        }
                      }}
                      role="button"
                      tabIndex={0}
                      aria-pressed={isSelected}
                      title={`Highlight CLR ${clrKey} in the 3D viewer`}
                      style={{
                        cursor: "pointer",
                        background: isSelected ? "#fff1ff" : "white",
                        outline: isSelected ? "2px solid #ff00ff" : "none",
                        outlineOffset: -2,
                      }}
                    >
                      <td style={styles.td}>
                        <strong>{r.rank}</strong>
                      </td>

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
                        <strong
                          style={{
                            color: isSelected ? "#b000b0" : "#111",
                            textDecoration: isSelected ? "underline" : "none",
                          }}
                        >
                          CLR {clrKey}
                        </strong>
                        {isSelected && (
                          <span style={styles.selectedBadge}>Selected</span>
                        )}
                      </td>

                      {modelOrder.map((modelName) => {
                        const data = r[modelName];
                        if (!data) return <td key={modelName} style={styles.tdMuted}>—</td>;

                        const score = Number(data.mean_score);
                        const scoreLabel = getScoreLabel(modelName, score);
                        const labelColor = getLabelColor(scoreLabel);

                        return (
                          <td key={modelName} style={styles.td}>
                            <div style={{ fontWeight: 700, color: labelColor }}>
                              {Number.isFinite(score) ? score.toFixed(4) : "—"}
                            </div>
                            <div
                              style={{
                                display: "inline-block",
                                marginTop: 3,
                                padding: "2px 7px",
                                borderRadius: 999,
                                border: `1px solid ${labelColor}`,
                                color: labelColor,
                                fontSize: 12,
                                fontWeight: 700,
                                lineHeight: 1.4,
                              }}
                            >
                              {scoreLabel}
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

          <div style={styles.interpretationSection}>
            <div style={styles.interpretationHeader}>
              <div>
                <h3 style={{ margin: 0 }}>Important residues and atoms</h3>
                <div style={styles.interpretationSubtitle}>
                  Select a CLR site, then choose an architecture to update the
                  structure highlighting.
                </div>
              </div>

              {selectedResult && (
                <div style={styles.modelSelector}>
                  {modelOrder.map((modelName) => {
                    const hasModelResult = Boolean(selectedResult?.[modelName]);
                    const hasInterpretation = Boolean(
                      selectedResult?.[modelName]?.interpretation
                    );
                    const active = selectedModel === modelName;
                    return (
                      <button
                        key={modelName}
                        type="button"
                        aria-pressed={active}
                        onMouseDown={(event) => event.stopPropagation()}
                        onClick={(event) => {
                          event.preventDefault();
                          event.stopPropagation();
                          setSelectedModel(modelName);
                        }}
                        title={
                          hasInterpretation
                            ? `Show ${modelName} interpretation`
                            : hasModelResult
                            ? `${modelName} score exists, but its interpretation payload is missing`
                            : `${modelName} model result is unavailable`
                        }
                        style={{
                          ...styles.modelButton,
                          ...(active ? styles.modelButtonActive : {}),
                          opacity: hasModelResult ? 1 : 0.65,
                        }}
                      >
                        {modelName}
                        {!hasInterpretation && (
                          <span aria-hidden="true" style={styles.modelWarningDot}>
                            !
                          </span>
                        )}
                      </button>
                    );
                  })}
                </div>
              )}
            </div>

            {!selectedResult && (
              <div style={styles.emptyInterpretation}>
                Click a CLR row above to view its model interpretation and exact
                PDB atom/residue IDs.
              </div>
            )}

            {selectedResult && selectedInterpretation && (
              <>
                <div style={styles.interpretationMeta}>
                  <strong>{selectedModel}</strong> · {selectedInterpretation.method}
                  <span style={styles.viewerLegendItem}>
                    <span
                      style={{ ...styles.legendSwatch, background: "#f59e0b" }}
                    />
                    important residue
                  </span>
                  <span style={styles.viewerLegendItem}>
                    <span
                      style={{ ...styles.legendSwatch, background: "#dc2626" }}
                    />
                    exact atom
                  </span>
                </div>

                {selectedInterpretation.warning && (
                  <div style={styles.interpretationWarning}>
                    {selectedInterpretation.warning}
                  </div>
                )}

                <div style={styles.interpretationGrid}>
                  <div style={styles.importanceCard}>
                    <h4 style={styles.importanceTitle}>Top residues</h4>
                    {Array.isArray(selectedInterpretation.top_residues) &&
                    selectedInterpretation.top_residues.length > 0 ? (
                      selectedInterpretation.top_residues.map((residue) => {
                        const importance = Number(residue.importance || 0);
                        const residueId = `${residue.chain_id || ""}${
                          residue.residue_number
                        }${residue.insertion_code || ""}`;
                        return (
                          <div
                            key={`${residueId}-${residue.residue_name}`}
                            style={styles.importanceItem}
                          >
                            <div style={styles.importanceLabel}>
                              <strong>{residue.residue_name}</strong> {residueId}
                            </div>
                            <div style={styles.importanceTrack}>
                              <div
                                style={{
                                  ...styles.importanceFillResidue,
                                  width: `${Math.max(
                                    2,
                                    Math.min(100, importance * 100)
                                  )}%`,
                                }}
                              />
                            </div>
                            <div style={styles.importanceValue}>
                              {(importance * 100).toFixed(0)}%
                            </div>
                          </div>
                        );
                      })
                    ) : (
                      <div style={styles.emptySmall}>
                        Exact residue IDs are unavailable for this site.
                      </div>
                    )}
                  </div>

                  <div style={styles.importanceCard}>
                    <h4 style={styles.importanceTitle}>Top atoms</h4>
                    {Array.isArray(selectedInterpretation.top_atoms) &&
                    selectedInterpretation.top_atoms.length > 0 ? (
                      selectedInterpretation.top_atoms.map((atom) => {
                        const importance = Number(atom.importance || 0);
                        const hasPdbId = Number.isFinite(Number(atom.atom_serial));
                        const residueId = `${atom.chain_id || ""}${
                          atom.residue_number ?? ""
                        }${atom.insertion_code || ""}`;
                        const label = hasPdbId
                          ? `${atom.residue_name} ${residueId} · ${
                              atom.atom_name || atom.atom_type
                            } · serial ${atom.atom_serial}`
                          : `graph node ${atom.node_index} · ${atom.atom_type}`;
                        return (
                          <div
                            key={`${atom.node_index}-${atom.atom_serial ?? "node"}`}
                            style={styles.importanceItem}
                          >
                            <div style={styles.importanceLabel}>{label}</div>
                            <div style={styles.importanceTrack}>
                              <div
                                style={{
                                  ...styles.importanceFillAtom,
                                  width: `${Math.max(
                                    2,
                                    Math.min(100, importance * 100)
                                  )}%`,
                                }}
                              />
                            </div>
                            <div style={styles.importanceValue}>
                              {(importance * 100).toFixed(0)}%
                            </div>
                          </div>
                        );
                      })
                    ) : (
                      <div style={styles.emptySmall}>
                        No positive atom-importance signal was found.
                      </div>
                    )}
                  </div>
                </div>
              </>
            )}

            {selectedResult && !selectedInterpretation && (
              <div style={styles.missingInterpretation}>
                <strong>{selectedModel} is selected.</strong>{" "}
                {selectedModelResult
                  ? "Its prediction score was returned, but its interpretation payload was not."
                  : "No prediction result was returned for this architecture."}
                {selectedModelResult && interpretationSchemaVersion === null && (
                  <>
                    {" "}The server appears to be using the older score-only API.
                    Restart the Python web service after deploying the updated
                    <code style={styles.inlineCode}> cholbind.py</code> and
                    <code style={styles.inlineCode}> application.py</code> files.
                  </>
                )}
                {selectedModelResult && interpretationSchemaVersion !== null && (
                  <>
                    {" "}Check the Python server log for an interpretation error or
                    a missing model checkpoint.
                  </>
                )}
              </div>
            )}
          </div>

          <div style={{ marginTop: 10, fontSize: 13, color: "#555" }}>
            Tip: click any evaluated CLR row to highlight and zoom to that exact CLR
            chain/residue ID. Orange marks important residues and red marks exact
            important atoms for the selected architecture.
          </div>
        </div>
      )}

      {/* ---------- Batch results ---------- */}
      {Array.isArray(batchResults) && (
        <div style={styles.resultsSection}>
          <h2>Batch Results</h2>
          <button style={styles.button} onClick={downloadBatchCsv}>
            Download CSV
          </button>

          <div style={{ marginTop: 10, fontSize: 13, color: "#555" }}>
            CSV contains one row per CLR per PDB file, including rank and the
            percentile-based label for every model score. Batch interpretations
            are aggregated by residue and atom type, without structure-specific IDs.
          </div>

          {batchInterpretation?.models && (
            <div style={styles.batchInterpretationGrid}>
              {modelOrder.map((modelName) => {
                const summary = batchInterpretation.models[modelName] || {};
                const residueTypes = summary.top_residue_types || [];
                const atomTypes = summary.top_atom_types || [];

                return (
                  <div key={modelName} style={styles.batchInterpretationCard}>
                    <h3 style={{ margin: 0 }}>{modelName}</h3>
                    <div style={styles.batchMethod}>
                      {summary.method || "Interpretation unavailable"}
                    </div>

                    <h4 style={styles.batchSubheading}>Important residue types</h4>
                    <div style={styles.typeChipWrap}>
                      {residueTypes.length > 0 ? (
                        residueTypes.map((record) => (
                          <span
                            key={record.residue_name}
                            style={styles.residueChip}
                            title={`${record.count} appearances among selected important atoms`}
                          >
                            {record.residue_name} · {Number(record.percent).toFixed(1)}%
                          </span>
                        ))
                      ) : (
                        <span style={styles.emptySmall}>No mapped residue types</span>
                      )}
                    </div>

                    <h4 style={styles.batchSubheading}>Important atom types</h4>
                    <div style={styles.typeChipWrap}>
                      {atomTypes.length > 0 ? (
                        atomTypes.map((record) => (
                          <span
                            key={record.atom_type}
                            style={styles.atomChip}
                            title={`${record.count} appearances among selected important atoms`}
                          >
                            {record.atom_type} · {Number(record.percent).toFixed(1)}%
                          </span>
                        ))
                      ) : (
                        <span style={styles.emptySmall}>No atom types available</span>
                      )}
                    </div>
                  </div>
                );
              })}
            </div>
          )}
        </div>
      )}
    </div>
  );
}

const styles = {
  body: {
    fontFamily: "sans-serif",
    maxWidth: "95vw",
    margin: "0 auto",
    padding: 20,
    lineHeight: 1.7,
    fontSize: 18,
  },
  h1: { color: "#333", fontSize: 34, marginBottom: 10 },
  uploadSection: {
    background: "#f9f9f9",
    padding: 22,
    borderRadius: 10,
    marginBottom: 12,
    border: "1px solid #ddd",
  },
  description: { fontSize: 18, color: "#333", marginBottom: 10 },
  uploadControls: {
    display: "flex",
    flexDirection: "column",
    alignItems: "center",
    gap: 14,
    marginTop: 10,
  },
  exampleRow: {
    display: "flex",
    alignItems: "center",
    justifyContent: "center",
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
  exampleHint: { color: "#555", fontSize: 14 },
  formRow: {
    display: "flex",
    alignItems: "center",
    justifyContent: "center",
    gap: 12,
    flexWrap: "wrap",
  },
  label: { fontWeight: 700, fontSize: 18 },
  fileInput: { fontSize: 16 },
  error: {
    color: "red",
    background: "#ffe6e6",
    padding: 10,
    borderRadius: 4,
    marginBottom: 12,
  },
  button: {
    padding: "10px 18px",
    cursor: "pointer",
    fontSize: 16,
    fontWeight: 700,
    borderRadius: 8,
    border: "1px solid #bbb",
    background: "white",
  },
  viewerNote: { marginTop: 12, fontSize: 16, color: "#444" },
  resultsSection: { marginTop: 20 },
  resultsHeader: {
    display: "flex",
    alignItems: "center",
    justifyContent: "space-between",
    gap: 12,
    flexWrap: "wrap",
  },
  selectedClrControls: {
    display: "flex",
    alignItems: "center",
    gap: 10,
    flexWrap: "wrap",
    fontSize: 14,
  },
  thresholdLegend: {
    marginTop: 12,
    padding: "10px 12px",
    border: "1px solid #ddd",
    borderRadius: 8,
    background: "#fafafa",
    color: "#444",
    fontSize: 13,
  },
  thresholdTitle: {
    marginBottom: 8,
    fontWeight: 700,
    color: "#333",
  },
  thresholdTable: {
    width: "100%",
    borderCollapse: "collapse",
    background: "white",
  },
  thresholdTh: {
    padding: "7px 9px",
    border: "1px solid #ddd",
    background: "#f3f3f3",
    textAlign: "left",
    whiteSpace: "nowrap",
  },
  thresholdTd: {
    padding: "7px 9px",
    border: "1px solid #ddd",
    whiteSpace: "nowrap",
  },
  selectedBadge: {
    display: "inline-block",
    marginLeft: 8,
    padding: "2px 7px",
    borderRadius: 999,
    background: "#ff00ff",
    color: "white",
    fontSize: 11,
    fontWeight: 700,
    verticalAlign: "middle",
  },
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
  interpretationSection: {
    marginTop: 18,
    padding: 16,
    border: "1px solid #ddd",
    borderRadius: 10,
    background: "#fafafa",
  },
  interpretationHeader: {
    display: "flex",
    alignItems: "center",
    justifyContent: "space-between",
    gap: 14,
    flexWrap: "wrap",
  },
  interpretationSubtitle: {
    marginTop: 3,
    color: "#555",
    fontSize: 14,
  },
  modelSelector: {
    display: "flex",
    gap: 7,
    flexWrap: "wrap",
    position: "relative",
    zIndex: 2,
    pointerEvents: "auto",
  },
  modelButton: {
    padding: "7px 13px",
    border: "1px solid #aaa",
    borderRadius: 999,
    background: "white",
    color: "#222",
    cursor: "pointer",
    fontSize: 14,
    fontWeight: 700,
    position: "relative",
    pointerEvents: "auto",
  },
  modelButtonActive: {
    background: "#111827",
    borderColor: "#111827",
    color: "white",
  },
  modelWarningDot: {
    display: "inline-flex",
    alignItems: "center",
    justifyContent: "center",
    width: 16,
    height: 16,
    marginLeft: 6,
    borderRadius: 999,
    background: "#b45309",
    color: "white",
    fontSize: 11,
    lineHeight: 1,
    verticalAlign: "middle",
  },
  missingInterpretation: {
    marginTop: 14,
    padding: "11px 13px",
    border: "1px solid #f0c36d",
    borderRadius: 8,
    background: "#fff8e5",
    color: "#6f4b00",
    fontSize: 14,
  },
  inlineCode: {
    padding: "1px 4px",
    borderRadius: 4,
    background: "rgba(255,255,255,0.75)",
  },
  emptyInterpretation: {
    marginTop: 14,
    padding: 16,
    borderRadius: 8,
    background: "white",
    color: "#555",
    fontSize: 15,
    textAlign: "center",
  },
  interpretationMeta: {
    display: "flex",
    alignItems: "center",
    gap: 14,
    flexWrap: "wrap",
    marginTop: 14,
    color: "#444",
    fontSize: 14,
  },
  viewerLegendItem: {
    display: "inline-flex",
    alignItems: "center",
    gap: 5,
  },
  legendSwatch: {
    display: "inline-block",
    width: 11,
    height: 11,
    borderRadius: 3,
    border: "1px solid rgba(0,0,0,0.25)",
  },
  interpretationWarning: {
    marginTop: 10,
    padding: "8px 10px",
    border: "1px solid #f0c36d",
    borderRadius: 7,
    background: "#fff8e5",
    color: "#6f4b00",
    fontSize: 13,
  },
  interpretationGrid: {
    display: "grid",
    gridTemplateColumns: "repeat(auto-fit, minmax(320px, 1fr))",
    gap: 14,
    marginTop: 12,
  },
  importanceCard: {
    padding: 13,
    border: "1px solid #ddd",
    borderRadius: 9,
    background: "white",
  },
  importanceTitle: {
    margin: "0 0 10px",
    color: "#222",
  },
  importanceItem: {
    display: "grid",
    gridTemplateColumns: "minmax(160px, 1.5fr) minmax(80px, 1fr) 42px",
    alignItems: "center",
    gap: 8,
    margin: "8px 0",
    fontSize: 13,
  },
  importanceLabel: {
    minWidth: 0,
    overflow: "hidden",
    textOverflow: "ellipsis",
    whiteSpace: "nowrap",
    color: "#333",
  },
  importanceTrack: {
    height: 8,
    overflow: "hidden",
    borderRadius: 999,
    background: "#e5e7eb",
  },
  importanceFillResidue: {
    height: "100%",
    borderRadius: 999,
    background: "#f59e0b",
  },
  importanceFillAtom: {
    height: "100%",
    borderRadius: 999,
    background: "#dc2626",
  },
  importanceValue: {
    color: "#555",
    textAlign: "right",
    fontVariantNumeric: "tabular-nums",
  },
  emptySmall: {
    color: "#777",
    fontSize: 13,
  },
  batchInterpretationGrid: {
    display: "grid",
    gridTemplateColumns: "repeat(auto-fit, minmax(280px, 1fr))",
    gap: 14,
    marginTop: 16,
  },
  batchInterpretationCard: {
    padding: 15,
    border: "1px solid #ddd",
    borderRadius: 10,
    background: "white",
  },
  batchMethod: {
    color: "#666",
    fontSize: 13,
    marginTop: 2,
  },
  batchSubheading: {
    margin: "13px 0 7px",
    fontSize: 15,
  },
  typeChipWrap: {
    display: "flex",
    flexWrap: "wrap",
    gap: 6,
  },
  residueChip: {
    display: "inline-block",
    padding: "4px 8px",
    border: "1px solid #d79520",
    borderRadius: 999,
    background: "#fff8e5",
    color: "#6f4b00",
    fontSize: 12,
    fontWeight: 700,
  },
  atomChip: {
    display: "inline-block",
    padding: "4px 8px",
    border: "1px solid #d66",
    borderRadius: 999,
    background: "#fff0f0",
    color: "#8f1d1d",
    fontSize: 12,
    fontWeight: 700,
  },

  publicationSection: {
  background: "white",
  padding: 20,
  marginBottom: 18,
  border: "1px solid #ddd",
  borderRadius: 10,
  textAlign: "center",
},

abstractFigureLink: {
  display: "block",
  width: "100%",
  textDecoration: "none",
},

abstractFigure: {
  display: "block",
  width: "100%",
  maxWidth: 1000,
  maxHeight: 500,
  objectFit: "contain",
  margin: "0 auto",
  borderRadius: 8,
},

publicationText: {
  maxWidth: 900,
  margin: "14px auto 8px",
  color: "#333",
  fontSize: 17,
  lineHeight: 1.6,
},

publicationLink: {
  display: "inline-block",
  marginTop: 6,
  color: "#065fd4",
  fontSize: 17,
  fontWeight: 700,
  textDecoration: "underline",
},
};