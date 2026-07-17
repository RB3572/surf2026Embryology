/* Sperm α/β analysis — 2-cell blastomere alpha/beta assignments across methods.
 * Recreates the mentor's grid (rows = methods, cols = embryos, cell = which original
 * blastomere is alpha: A = labels 1+3, B = labels 2+4). Rows can be flipped (the α/β
 * direction is arbitrary per method); Auto-align flips rows to maximise column agreement.
 * Data: build_alphabeta.py -> data/alphabeta.json. Front-end only reads it.
 */
(() => {
  const $ = (s) => document.querySelector(s);
  const el = (tag, cls) => { const e = document.createElement(tag); if (cls) e.className = cls; return e; };
  const STAGES = ["early", "late"];
  const state = { data: null, flips: { early: [], late: [] } };

  const flipCall = (c) => (c === "A" ? "B" : c === "B" ? "A" : c);            // T / NA unchanged
  const shown = (call, flipped) => (flipped ? flipCall(call) : call);

  function consensus(calls) {
    let a = 0, b = 0;
    for (const c of calls) { if (c === "A") a++; else if (c === "B") b++; }
    if (a === 0 && b === 0) return "—";
    return a > b ? "A" : b > a ? "B" : "T";
  }

  // ---------- render ----------
  function renderLegend() {
    const L = state.data.legend, host = $("#ab-legend");
    const items = [["A", L.A], ["B", L.B], ["T", L.T], ["NA", L.NA]];
    host.innerHTML = items.map(([k, txt]) =>
      `<span class="ab-key"><span class="ab-sw c-${k}">${k === "NA" ? "" : k}</span>${txt}</span>`).join("");
  }

  function tooltip(r, m, disp) {
    const raw = r.raw;
    let s = `${r.id}\n${m.label}: ${disp === "NA" ? "unavailable" : disp === "T" ? "tie" : "α = blastomere " + disp}`;
    if (m.id === "total") s += `\nblastomere A (1+3) = ${raw.nA.toLocaleString()} transcripts\nblastomere B (2+4) = ${raw.nB.toLocaleString()}`;
    if (m.id === "sperm") s += raw.spermSeg != null ? `\nsperm in segment ${raw.spermSeg}` : `\nno labelled sperm in this embryo`;
    return s;
  }

  function renderGrid(stage) {
    const grid = $(`#grid-${stage}`), rows = state.data.stages[stage], methods = state.data.methods, flips = state.flips[stage];
    grid.style.setProperty("--cols", rows.length);
    grid.innerHTML = "";
    const frag = document.createDocumentFragment();

    methods.forEach((m, mi) => {
      const lab = el("div", "ab-rowlabel");
      const btn = el("button", "ab-flip" + (flips[mi] ? " on" : ""));
      btn.textContent = "⇄"; btn.title = "Flip α/β for this row";
      btn.addEventListener("click", () => { flips[mi] = !flips[mi]; renderGrid(stage); });
      const name = el("span", "ab-rowname"); name.textContent = m.label;
      lab.appendChild(btn); lab.appendChild(name);
      frag.appendChild(lab);
      rows.forEach((r) => {
        const disp = shown(r.calls[m.id], flips[mi]);
        const cell = el("div", `ab-cell c-${disp}`);
        cell.textContent = disp;
        cell.title = tooltip(r, m, disp);
        frag.appendChild(cell);
      });
    });

    // consensus (derived) row
    const clab = el("div", "ab-rowlabel ab-consensus-lab"); clab.innerHTML = "<span class='ab-rowname'>Consensus</span>";
    frag.appendChild(clab);
    rows.forEach((r) => {
      const calls = methods.map((m, mi) => shown(r.calls[m.id], flips[mi]));
      const c = consensus(calls);
      const cell = el("div", `ab-cell ab-cons c-${c === "—" ? "NA" : c}`);
      cell.textContent = c; cell.title = `${r.id}\nconsensus α = ${c}`;
      frag.appendChild(cell);
    });

    // embryo id labels (rotated)
    frag.appendChild(el("div", "ab-corner"));
    rows.forEach((r) => {
      const cl = el("div", "ab-collabel"); const sp = document.createElement("span");
      sp.textContent = r.id; cl.appendChild(sp); frag.appendChild(cl);
    });
    grid.appendChild(frag);
  }

  function renderAll() { STAGES.forEach(renderGrid); }

  // ---------- flip operations ----------
  function autoAlignStage(stage) {
    const rows = state.data.stages[stage], methods = state.data.methods;
    const ne = rows.length, nm = methods.length;
    const M = methods.map((m) => rows.map((r) => { const c = r.calls[m.id]; return c === "A" ? 1 : c === "B" ? -1 : 0; }));
    const flip = state.flips[stage].slice();
    for (let iter = 0; iter < 100; iter++) {
      let changed = false;
      for (let r = 0; r < nm; r++) {
        let same = 0, opp = 0;
        for (let c = 0; c < ne; c++) {
          if (M[r][c] === 0) continue;
          let cons = 0;
          for (let r2 = 0; r2 < nm; r2++) { if (r2 === r) continue; cons += (flip[r2] ? -1 : 1) * M[r2][c]; }
          if (cons === 0) continue;
          if (Math.sign(M[r][c]) === Math.sign(cons)) same++; else opp++;
        }
        const want = opp > same;
        if (want !== flip[r]) { flip[r] = want; changed = true; }
      }
      if (!changed) break;
    }
    state.flips[stage] = flip;
  }
  function autoAlign() { STAGES.forEach(autoAlignStage); renderAll(); }
  function flipAll() { STAGES.forEach((s) => { state.flips[s] = state.flips[s].map((f) => !f); }); renderAll(); }
  function resetFlips() { STAGES.forEach((s) => { state.flips[s] = state.data.methods.map(() => false); }); renderAll(); }

  // ---------- CSV ----------
  function exportCSV() {
    const methods = state.data.methods;
    const head = ["stage", "embryo", ...methods.map((m) => m.label.replace(/,/g, "")), "consensus_alpha", "nA_labels_1_3", "nB_labels_2_4", "sperm_segment"];
    const lines = [head.join(",")];
    for (const stage of STAGES) {
      const rows = state.data.stages[stage], flips = state.flips[stage];
      rows.forEach((r) => {
        const calls = methods.map((m, mi) => shown(r.calls[m.id], flips[mi]));
        lines.push([stage, r.id, ...calls, consensus(calls), r.raw.nA, r.raw.nB, r.raw.spermSeg == null ? "" : r.raw.spermSeg].join(","));
      });
    }
    const blob = new Blob([lines.join("\n")], { type: "text/csv" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob); a.download = "sperm_alpha_beta_2cell.csv";
    document.body.appendChild(a); a.click(); a.remove(); URL.revokeObjectURL(a.href);
  }

  // ---------- boot ----------
  (async function init() {
    try {
      state.data = await (await fetch("data/alphabeta.json")).json();
    } catch (err) {
      $("#ab-loading").textContent = "Failed to load data: " + (err.message || err); return;
    }
    $("#ab-loading").hidden = true;
    STAGES.forEach((s) => { state.flips[s] = state.data.methods.map(() => false); });
    const nE = state.data.stages.early.length, nL = state.data.stages.late.length;
    const spE = state.data.stages.early.filter((r) => r.calls.sperm !== "NA").length;
    const spL = state.data.stages.late.filter((r) => r.calls.sperm !== "NA").length;
    $("#embryo-count").textContent = `${nE + nL} embryos · ${spE + spL} with sperm`;
    $("#sub-early").textContent = `${nE} embryos · ${spE} with a labelled sperm`;
    $("#sub-late").textContent = `${nL} embryos · ${spL} with a labelled sperm`;
    renderLegend(); renderAll();
    $("#ab-auto").addEventListener("click", autoAlign);
    $("#ab-flipall").addEventListener("click", flipAll);
    $("#ab-reset").addEventListener("click", resetFlips);
    $("#ab-csv").addEventListener("click", exportCSV);
  })();
})();
