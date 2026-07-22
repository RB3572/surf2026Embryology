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

  // ---- project access matrix ----
  async function loadAccess() {
    const el = $("#ad-access"), note = $("#ad-access-note");
    let d;
    try { d = await (await fetch("/api/access?matrix=1", { credentials: "same-origin", cache: "no-store" })).json(); }
    catch (e) { note.textContent = "failed to load: " + e.message; return; }
    if (!d.ok) { el.innerHTML = ""; note.textContent = "access config unavailable — " + (d.err || "database error"); return; }
    const projs = d.projects || [];
    const head = `<thead><tr><th class="ad-first">Login</th>${projs.map((p) => `<th title="${esc(p.label)}">${esc(p.label)}</th>`).join("")}</tr></thead>`;
    const body = (d.users || []).map((u) => {
      const isAdmin = u.role === "admin", all = u.projects == null;
      const set = new Set(u.projects || []);
      const cells = projs.map((p) => {
        const on = isAdmin || all || set.has(p.key);
        return `<td class="ad-acc"><input type="checkbox" data-usr="${esc(u.usr)}" data-proj="${p.key}"${on ? " checked" : ""}${isAdmin ? " disabled" : ""}></td>`;
      }).join("");
      return `<tr><td class="ad-first"><b>${esc(u.usr)}</b>${isAdmin ? ' <span class="ad-dim">(admin · all)</span>' : all ? ' <span class="ad-dim">(all)</span>' : ""}</td>${cells}</tr>`;
    }).join("");
    el.innerHTML = head + `<tbody>${body}</tbody>`;
    note.textContent = "Uncheck a project to hide it from that login (they are redirected to the landing). Admins always have full access.";
    el.querySelectorAll("input[type=checkbox]").forEach((cb) => cb.addEventListener("change", onAccessToggle));
  }
  async function onAccessToggle(e) {
    const usr = e.target.dataset.usr;
    const checked = [...$("#ad-access").querySelectorAll(`input[data-usr="${CSS.escape(usr)}"]`)].filter((c) => c.checked).map((c) => c.dataset.proj);
    $("#ad-access-note").textContent = "saving…";
    try {
      const d = await (await fetch("/api/access", { method: "POST", credentials: "same-origin",
        headers: { "content-type": "application/json" }, body: JSON.stringify({ usr, projects: checked }) })).json();
      $("#ad-access-note").textContent = d.ok ? `Saved — ${usr} now has ${checked.length} project(s).` : "save failed: " + (d.err || "unknown");
    } catch (err) { $("#ad-access-note").textContent = "save failed: " + err.message; }
  }

  // ---- logins & passwords ----
  async function loadUsers() {
    const el = $("#ad-logins"), note = $("#ad-adduser-note");
    let d;
    try { d = await (await fetch("/api/users", { credentials: "same-origin", cache: "no-store" })).json(); }
    catch (e) { el.innerHTML = ""; note.textContent = "failed to load: " + e.message; return; }
    if (!d.ok) { el.innerHTML = ""; note.textContent = "logins unavailable — " + (d.err || "error"); return; }
    const rows = (d.accounts || []).map((a) => {
      const role = a.role === "admin" ? `<span class="ad-role-admin">admin</span>` : `<span class="ad-dim">user</span>`;
      const last = a.kind === "added"
        ? `<button type="button" class="ad-rm" data-usr="${esc(a.usr)}" title="Remove this login">Remove</button>`
        : `<span class="ad-dim">built-in</span>`;
      return [`<b>${esc(a.usr)}</b>`, role, `<code class="ad-pw">${esc(a.pw)}</code>`, last];
    });
    table(el, ["Login", "Role", "Password", ""], rows);
    el.querySelectorAll(".ad-rm").forEach((b) => b.addEventListener("click", () => removeUser(b.dataset.usr)));
    if (d.hasDb === false && !note.textContent)
      note.textContent = "No database configured — added logins won't persist (built-in logins still work).";
  }
  async function removeUser(usr) {
    if (!window.confirm(`Remove the login "${usr}"? They will no longer be able to sign in.`)) return;
    const note = $("#ad-adduser-note"); note.textContent = "removing…";
    try {
      const d = await (await fetch("/api/users?usr=" + encodeURIComponent(usr), { method: "DELETE", credentials: "same-origin" })).json();
      note.textContent = d.ok ? `Removed “${usr}”.` : "remove failed: " + (d.err || "unknown");
      if (d.ok) { loadUsers(); loadAccess(); }
    } catch (e) { note.textContent = "remove failed: " + e.message; }
  }
  $("#ad-adduser").addEventListener("submit", async (e) => {
    e.preventDefault();
    const usr = $("#ad-nu-name").value.trim(), pw = $("#ad-nu-pw").value, role = $("#ad-nu-role").value;
    const note = $("#ad-adduser-note");
    if (!usr || !pw) { note.textContent = "enter a login name and a password"; return; }
    note.textContent = "adding…";
    try {
      const d = await (await fetch("/api/users", { method: "POST", credentials: "same-origin",
        headers: { "content-type": "application/json" }, body: JSON.stringify({ usr, pw, role }) })).json();
      if (d.ok) {
        note.textContent = `Added “${usr}” — they can sign in with that password now.`;
        $("#ad-nu-name").value = ""; $("#ad-nu-pw").value = "";
        loadUsers(); loadAccess();
      } else note.textContent = "couldn't add: " + (d.err || "unknown");
    } catch (err) { note.textContent = "couldn't add: " + err.message; }
  });

  loadUsers();
  loadAccess();

  $("#ad-user-filter").addEventListener("change", (e) => setFilter(e.target.value));
  $("#ad-users").addEventListener("click", (e) => {
    const tr = e.target.closest("tbody tr"); if (!tr) return;
    const name = (tr.querySelector(".ad-first") || {}).textContent;
    if (name && name !== "—") setFilter(name === filterUser ? "" : name);   // click the active row again to clear
  });
  $("#ad-refresh").addEventListener("click", load);
  load();
})();
