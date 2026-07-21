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

  let filterUser = "";   // "" = all users; otherwise a single login name

  function populateFilter(users) {
    const sel = $("#ad-user-filter");
    const names = users.map((u) => u.usr).filter(Boolean);
    if (filterUser && !names.includes(filterUser)) names.unshift(filterUser);   // keep an out-of-range choice
    sel.innerHTML = `<option value="">All users</option>` +
      names.map((n) => `<option value="${esc(n)}"${n === filterUser ? " selected" : ""}>${esc(n)}</option>`).join("");
    sel.value = filterUser;
  }
  function renderFilterBar() {
    const bar = $("#ad-filterbar");
    if (!filterUser) { bar.hidden = true; bar.innerHTML = ""; return; }
    bar.hidden = false;
    bar.innerHTML = `<span>Filtered to <b>${esc(filterUser)}</b> — projects, genes, downloads and activity below are for this login only.</span>` +
      `<button type="button" class="ad-clear" id="ad-clear">Show all users ✕</button>`;
    $("#ad-clear").addEventListener("click", () => setFilter(""));
  }
  function setFilter(u) { if (u === filterUser) return; filterUser = u; $("#ad-user-filter").value = u; load(); }

  async function load() {
    $("#ad-status").textContent = "loading…";
    let d;
    try {
      const url = "/api/admin" + (filterUser ? "?user=" + encodeURIComponent(filterUser) : "");
      const r = await fetch(url, { credentials: "same-origin", cache: "no-store" });
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

    populateFilter(d.users || []);
    renderFilterBar();

    const t = d.totals || {};
    $("#ad-kpis").innerHTML = [
      ["People", filterUser ? 1 : (d.users || []).length],
      [filterUser ? "Their events" : "Total events", t.events || 0],
      ["Downloads", (d.downloads || []).length],
      ["Genes tracked", (d.genes || []).length],
    ].map(([k, v]) => `<div class="ad-kpi"><div class="ad-kpi-n">${esc(v)}</div><div class="ad-kpi-k">${esc(k)}</div></div>`).join("");

    table($("#ad-users"), ["Person", "Views", "Downloads", "Projects", "Active days", "Events", "Last seen"],
      (d.users || []).map((u) => [
        `<b>${esc(u.usr || "—")}</b>`, u.views || 0, u.downloads || 0, u.projects || 0,
        u.active_days || 0, u.events || 0, `${esc(fmtTs(u.last_seen))} <span class="ad-dim">(${ago(u.last_seen)})</span>`]));
    $("#ad-users").classList.add("ad-clickable");   // click a person to filter to them
    $("#ad-users").querySelectorAll("tbody tr").forEach((tr) => {
      const n = (tr.querySelector(".ad-first") || {}).textContent;
      tr.classList.toggle("ad-active", !!filterUser && n === filterUser);
    });

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

  $("#ad-user-filter").addEventListener("change", (e) => setFilter(e.target.value));
  $("#ad-users").addEventListener("click", (e) => {
    const tr = e.target.closest("tbody tr"); if (!tr) return;
    const name = (tr.querySelector(".ad-first") || {}).textContent;
    if (name && name !== "—") setFilter(name === filterUser ? "" : name);   // click the active row again to clear
  });
  $("#ad-refresh").addEventListener("click", load);
  load();
})();
