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

  // "" = everyone; otherwise a person KEY (a Supabase user id, or a display name for rows
  // recorded before the SSO migration).
  let filterUser = "";
  let filterLabel = "";

  function populateFilter(users) {
    const sel = $("#ad-user-filter");
    const opts = users.map((u) => ({ key: u.person, label: u.usr || u.email || u.person }))
      .filter((o) => o.key);
    if (filterUser && !opts.some((o) => o.key === filterUser)) {
      opts.unshift({ key: filterUser, label: filterLabel || filterUser });   // keep an out-of-range choice
    }
    sel.innerHTML = `<option value="">All users</option>` +
      opts.map((o) => `<option value="${esc(o.key)}"${o.key === filterUser ? " selected" : ""}>${esc(o.label)}</option>`).join("");
    sel.value = filterUser;
  }
  function renderFilterBar() {
    const bar = $("#ad-filterbar");
    if (!filterUser) { bar.hidden = true; bar.innerHTML = ""; return; }
    bar.hidden = false;
    bar.innerHTML = `<span>Filtered to <b>${esc(filterLabel || filterUser)}</b> — projects, genes, downloads and activity below are for this person only.</span>` +
      `<button type="button" class="ad-clear" id="ad-clear">Show all users ✕</button>`;
    $("#ad-clear").addEventListener("click", () => setFilter(""));
  }
  function setFilter(u, label) {
    if (u === filterUser) return;
    filterUser = u; filterLabel = label || "";
    $("#ad-user-filter").value = u;
    load();
  }

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

    const people = d.users || [];
    table($("#ad-users"), ["Person", "Views", "Downloads", "Projects", "Active days", "Events", "Last seen"],
      people.map((u) => [
        `<b>${esc(u.usr || u.email || "—")}</b>` +
          (u.email && u.usr ? ` <span class="ad-dim">${esc(u.email)}</span>` : ""),
        u.views || 0, u.downloads || 0, u.projects || 0,
        u.active_days || 0, u.events || 0, `${esc(fmtTs(u.last_seen))} <span class="ad-dim">(${ago(u.last_seen)})</span>`]));
    $("#ad-users").classList.add("ad-clickable");   // click a person to filter to them
    // Carry the person KEY on the row — display names are not unique and are no longer the key.
    $("#ad-users").querySelectorAll("tbody tr").forEach((tr, i) => {
      const u = people[i];
      if (!u) return;
      tr.dataset.person = u.person || "";
      tr.dataset.label = u.usr || u.email || u.person || "";
      tr.classList.toggle("ad-active", !!filterUser && u.person === filterUser);
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

  // ---- project access matrix ----
  async function loadAccess() {
    const el = $("#ad-access"), note = $("#ad-access-note");
    let d;
    try { d = await (await fetch("/api/access?matrix=1", { credentials: "same-origin", cache: "no-store" })).json(); }
    catch (e) { note.textContent = "failed to load: " + e.message; return; }
    if (!d.ok) { el.innerHTML = ""; note.textContent = "access config unavailable — " + (d.err || "database error"); return; }
    const projs = d.projects || [];
    const head = `<thead><tr><th class="ad-first">Member</th>${projs.map((p) => `<th title="${esc(p.label)}">${esc(p.label)}</th>`).join("")}</tr></thead>`;
    const body = (d.users || []).map((u) => {
      const all = u.projects == null;
      const set = new Set(u.projects || []);
      const cells = projs.map((p) => {
        const on = all || set.has(p.key);
        return `<td class="ad-acc"><input type="checkbox" data-key="${esc(u.key)}" data-proj="${p.key}"${on ? " checked" : ""}></td>`;
      }).join("");
      const who = esc(u.name || u.email || u.user_id || u.key);
      const tag = u.pending
        ? ' <span class="ad-dim">(pending first sign-in)</span>'
        : all ? ' <span class="ad-dim">(all)</span>' : "";
      const sub = u.name && u.email ? `<br><span class="ad-dim">${esc(u.email)}</span>` : "";
      return `<tr><td class="ad-first"><b>${who}</b>${tag}${sub}</td>${cells}</tr>`;
    }).join("");
    el.innerHTML = head + `<tbody>${body}</tbody>`;
    note.textContent = "Uncheck a project to hide it from that member (they are redirected to the landing). " +
      "Ticking every box removes the restriction. Admins and the PI always have full access.";
    el.querySelectorAll("input[type=checkbox]").forEach((cb) => cb.addEventListener("change", onAccessToggle));
  }
  async function onAccessToggle(e) {
    const key = e.target.dataset.key;
    const boxes = [...$("#ad-access").querySelectorAll(`input[data-key="${CSS.escape(key)}"]`)];
    const checked = boxes.filter((c) => c.checked).map((c) => c.dataset.proj);
    const note = $("#ad-access-note");
    note.textContent = "saving…";
    try {
      // All boxes ticked = no restriction at all, so drop the row entirely.
      const clearAll = checked.length === boxes.length;
      const d = clearAll
        ? await (await fetch("/api/access?key=" + encodeURIComponent(key),
            { method: "DELETE", credentials: "same-origin" })).json()
        : await (await fetch("/api/access", { method: "POST", credentials: "same-origin",
            headers: { "content-type": "application/json" },
            body: JSON.stringify({ key, projects: checked }) })).json();
      note.textContent = d.ok
        ? (clearAll ? "Saved — full access restored." : `Saved — ${checked.length} project(s).`)
        : "save failed: " + (d.err || "unknown");
    } catch (err) { note.textContent = "save failed: " + err.message; }
  }

  // Pre-restrict somebody by email, before they have ever signed in. The row is rebound to
  // their real Supabase user id on first sign-in, so the limit applies from their first visit.
  $("#ad-preseed").addEventListener("submit", async (e) => {
    e.preventDefault();
    const email = $("#ad-ps-email").value.trim();
    const note = $("#ad-preseed-note");
    if (!email) { note.textContent = "enter an email"; return; }
    note.textContent = "adding…";
    try {
      // Start with no projects; the admin then ticks what they should see.
      const d = await (await fetch("/api/access", { method: "POST", credentials: "same-origin",
        headers: { "content-type": "application/json" },
        body: JSON.stringify({ email, projects: [] }) })).json();
      if (d.ok) {
        note.textContent = `Added ${email} — now tick the projects they should see.`;
        $("#ad-ps-email").value = "";
        loadAccess();
      } else note.textContent = "couldn't add: " + (d.err || "unknown");
    } catch (err) { note.textContent = "couldn't add: " + err.message; }
  });

  // NOTE: there is no login/password management any more. Who may sign in is decided entirely
  // by Lab Logger lab membership (checked against Supabase at sign-in), so this console only
  // narrows what an already-authorised member can open.

  loadAccess();

  $("#ad-user-filter").addEventListener("change", (e) => {
    const sel = e.target;
    setFilter(sel.value, sel.options[sel.selectedIndex] ? sel.options[sel.selectedIndex].text : "");
  });
  $("#ad-users").addEventListener("click", (e) => {
    const tr = e.target.closest("tbody tr"); if (!tr) return;
    const key = tr.dataset.person;
    if (key) setFilter(key === filterUser ? "" : key, tr.dataset.label);   // click again to clear
  });
  $("#ad-refresh").addEventListener("click", load);
  load();
})();
