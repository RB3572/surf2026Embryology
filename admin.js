/* Admin usage console. Reads /api/admin (admin-gated) and renders the tables. Read-only. */
(() => {
  const $ = (s) => document.querySelector(s);
  const esc = (s) => String(s == null ? "" : s).replace(/[&<>"]/g, (c) =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const fmtTs = (t) => { try { return new Date(t).toLocaleString(); } catch { return t; } };
  const ago = (t) => {
    const s = (Date.now() - new Date(t).getTime()) / 1000;
    if (isNaN(s)) return "";
    if (s < 60) return Math.round(s) + "s ago";
    if (s < 3600) return Math.round(s / 60) + "m ago";
    if (s < 86400) return Math.round(s / 3600) + "h ago";
    return Math.round(s / 86400) + "d ago";
  };
  const table = (el, head, rows) => {
    el.innerHTML = `<thead><tr>${head.map((h) => `<th>${esc(h)}</th>`).join("")}</tr></thead>` +
      `<tbody>${rows.length ? rows.map((r) => `<tr>${r.map((c, i) =>
        `<td${i === 0 ? " class='ad-first'" : ""}>${c}</td>`).join("")}</tr>`).join("")
        : `<tr><td class="ad-empty" colspan="${head.length}">no data yet</td></tr>`}</tbody>`;
  };

  async function load() {
    $("#ad-status").textContent = "loading…";
    let d;
    try {
      const r = await fetch("/api/admin", { credentials: "same-origin", cache: "no-store" });
      d = await r.json();
    } catch (e) { $("#ad-status").textContent = "request failed: " + e.message; return; }

    if (d.ok === false) {
      $("#ad-status").className = "ad-status ad-err";
      $("#ad-status").textContent = "database error — " + (d.err || "unknown");
    } else if (d.empty || !(d.totals && d.totals.events)) {
      $("#ad-status").className = "ad-status";
      $("#ad-status").textContent = "connected · no activity recorded yet";
    } else {
      $("#ad-status").className = "ad-status ad-ok";
      $("#ad-status").textContent = "connected · live";
    }

    const t = d.totals || {};
    $("#ad-kpis").innerHTML = [
      ["People", (d.users || []).length],
      ["Total events", t.events || 0],
      ["Downloads", (d.downloads || []).length],
      ["Genes tracked", (d.genes || []).length],
    ].map(([k, v]) => `<div class="ad-kpi"><div class="ad-kpi-n">${esc(v)}</div><div class="ad-kpi-k">${esc(k)}</div></div>`).join("");

    table($("#ad-users"), ["Person", "Views", "Downloads", "Projects", "Active days", "Events", "Last seen"],
      (d.users || []).map((u) => [
        `<b>${esc(u.usr || "—")}</b>`, u.views || 0, u.downloads || 0, u.projects || 0,
        u.active_days || 0, u.events || 0, `${esc(fmtTs(u.last_seen))} <span class="ad-dim">(${ago(u.last_seen)})</span>`]));

    table($("#ad-projects"), ["Project", "Views", "Events", "People"],
      (d.projects || []).map((p) => [esc(p.project || "—"), p.views || 0, p.events || 0, p.users || 0]));

    table($("#ad-genes"), ["Gene", "Views", "People"],
      (d.genes || []).map((g) => [`<b>${esc(g.gene)}</b>`, g.n || 0, g.users || 0]));

    table($("#ad-downloads"), ["When", "Person", "Project", "What"],
      (d.downloads || []).map((x) => [esc(fmtTs(x.ts)), esc(x.usr), esc(x.project),
        esc((x.detail && (x.detail.control || "")) + (x.detail && x.detail.gene ? " · " + x.detail.gene : ""))]));

    table($("#ad-recent"), ["When", "Person", "Project", "Action", "Detail"],
      (d.recent || []).map((e) => {
        const det = e.detail || {};
        const s = det.gene ? "gene " + det.gene : det.control ? det.control : det.title ? det.title : (e.path || "");
        return [`${esc(fmtTs(e.ts))} <span class="ad-dim">(${ago(e.ts)})</span>`, esc(e.usr),
          esc(e.project), `<span class="ad-act ad-act-${esc(e.action)}">${esc(e.action)}</span>`, esc(s)];
      }));
  }

  $("#ad-refresh").addEventListener("click", load);
  load();
})();
