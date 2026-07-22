// Password-gate the whole site at the edge (Vercel Routing/Edge Middleware).
//
// Runs BEFORE any file is served, so nothing — pages, scripts, or the /data files — is
// reachable without a valid password. Framework-agnostic: plain Web APIs, no imports.
// Returning nothing continues to the static file; returning a Response short-circuits.
//
// Each password maps to an identity. `token` is the opaque cookie value (server-side only,
// never shipped to the browser, unguessable) so a session can be attributed to a person for
// analytics. `role: "admin"` unlocks the admin console; every /admin* path is 404 (not 403)
// for everyone else, so its existence is never even revealed.

import { neon } from "@neondatabase/serverless";

const COOKIE = "surf_gate";
const MAX_AGE = 60 * 60 * 24 * 30; // remember a login for 30 days

const ACCOUNTS = [
  { pw: "Sonichedgehog1", token: "s1_kathytam_b7f3a92c8d1e4056a19d", user: "Kathy Tam", role: "user" },
  { pw: "AdminPass",      token: "adm_owner_5c1e9a2f7b3d8064c2a7f1", user: "Admin",     role: "admin" },
];
const BY_TOKEN = new Map(ACCOUNTS.map((a) => [a.token, a]));

// ---- per-user project access (managed from the admin console; see api/access.mjs) ----
// A non-admin user restricted to a project list is redirected to the landing for any other
// project page. The list is read from Neon and cached; a DB failure fails OPEN (allow) so a
// database hiccup can never lock authenticated users out of the whole site.
const CONN = process.env.DATABASE_URL || process.env.POSTGRES_URL ||
             process.env.DATABASE_URL_UNPOOLED || process.env.POSTGRES_URL_NON_POOLING;
const accessSql = CONN ? neon(CONN) : null;
const PROJECT_OF = {
  "pronuclei.html": "pronuclei", "extpt.html": "extpt", "segments.html": "segments",
  "zygote-planes.html": "zygote-planes", "sperm-division.html": "sperm-division",
  "pronuclei-assignments.html": "pronuclei-assignments", "sperm-pca.html": "sperm-pca",
  "diffusion.html": "diffusion",
  "axes.html": "axes", "alphabeta.html": "alphabeta",
};
let _accMap = null, _accAt = 0;
async function allowedProjects(user) {          // → array of allowed keys, or null = all projects
  const now = Date.now();
  if (accessSql && (!_accMap || now - _accAt > 60000)) {
    try {
      const rows = await accessSql`SELECT usr, projects FROM access`;
      const m = {}; rows.forEach((r) => (m[r.usr] = r.projects));
      _accMap = m; _accAt = now;
    } catch (_) { if (!_accMap) _accMap = {}; }  // keep last-known; empty = everyone sees all
  } else if (!accessSql && !_accMap) { _accMap = {}; }
  const v = _accMap && _accMap[user];
  return Array.isArray(v) ? v : null;
}

export const config = {
  matcher: "/((?!_vercel).*)", // gate every path except Vercel's own internals
};

function cookieValue(cookies, name) {
  for (const c of cookies.split(/;\s*/)) {
    const i = c.indexOf("=");
    if (i > 0 && c.slice(0, i) === name) return c.slice(i + 1);
  }
  return null;
}

export default async function middleware(request) {
  const url = new URL(request.url);
  const cookies = request.headers.get("cookie") || "";
  const account = BY_TOKEN.get(cookieValue(cookies, COOKIE));

  // Log out: clear the session cookies and bounce to the gate.
  if (url.pathname === "/logout") {
    const headers = new Headers({ location: "/", "cache-control": "no-store" });
    headers.append("set-cookie", `${COOKIE}=; Path=/; Max-Age=0; HttpOnly; Secure; SameSite=Lax`);
    headers.append("set-cookie", `surf_admin=; Path=/; Max-Age=0; Secure; SameSite=Lax`);
    return new Response(null, { status: 303, headers });
  }

  // Admin-only surface. Anything under /admin (page, assets, data) is 404 unless this is an
  // admin session — indistinguishable from "does not exist", so non-admins never see it.
  const p = url.pathname;
  const adminOnly = p === "/admin" || p.startsWith("/admin/") || p.startsWith("/admin.") ||
                    p.startsWith("/api/admin");
  if (adminOnly && (!account || account.role !== "admin")) {
    return new Response("Not Found", {
      status: 404,
      headers: { "content-type": "text/plain; charset=utf-8", "cache-control": "no-store" },
    });
  }

  // Authenticated → project-access check, then serve. (The admin console card is revealed on the
  // landing client-side from the gated /admin.card.html fragment — see index.html — because an
  // edge self-fetch to inject it deadlocks on a custom domain. The fragment + page stay 404 for
  // everyone else, so the console is still invisible to non-admins.)
  if (account) {
    if (account.role !== "admin") {
      const proj = PROJECT_OF[p.replace(/^\//, "")];
      if (proj) {
        try {
          const allowed = await allowedProjects(account.user);
          if (allowed && !allowed.includes(proj)) {
            const to = new URL("/", url); to.searchParams.set("denied", proj);
            return Response.redirect(to.toString(), 303);
          }
        } catch (_) { /* fail-open: serve */ }
      }
    }
    return;
  }

  // login form submitted
  if (request.method === "POST") {
    let pwd = "";
    try { pwd = String((await request.formData()).get("password") || ""); } catch (_) { /* malformed */ }
    const match = ACCOUNTS.find((a) => a.pw === pwd);
    if (match) {
      const headers = new Headers({ location: url.pathname + url.search, "cache-control": "no-store" });
      headers.append("set-cookie", `${COOKIE}=${match.token}; Path=/; Max-Age=${MAX_AGE}; HttpOnly; Secure; SameSite=Lax`);
      // readable (non-HttpOnly) hint so the landing knows to pull in the admin card; useless to
      // forge — the card fragment + console are still server-gated to the admin token.
      if (match.role === "admin") headers.append("set-cookie", `surf_admin=1; Path=/; Max-Age=${MAX_AGE}; Secure; SameSite=Lax`);
      return new Response(null, { status: 303, headers });
    }
    return gate(true); // wrong password
  }

  return gate(false);
}

function gate(wrong) {
  return new Response(gateHtml(wrong), {
    status: 401,
    headers: { "content-type": "text/html; charset=utf-8", "cache-control": "no-store" },
  });
}

// Minimal black-and-white gate: a password box and a submit button, nothing else.
function gateHtml(wrong) {
  return `<!doctype html><html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Password</title>
<style>
  :root { color-scheme: light; }
  * { box-sizing: border-box; }
  html, body { height: 100%; margin: 0; }
  body { background: #fff; color: #000; display: grid; place-items: center; padding: 24px;
    font: 15px/1.4 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }
  form { display: flex; flex-direction: column; gap: 12px; width: 240px; }
  input, button { font-size: 15px; padding: 11px 12px; border: 1px solid #000; border-radius: 0;
    background: #fff; color: #000; outline: none; }
  input:focus { border-width: 2px; padding: 10px 11px; }
  button { cursor: pointer; font-weight: 600; }
  button:hover, button:focus { background: #000; color: #fff; }
  form.err input { animation: shake .28s; }
  @keyframes shake { 0%,100%{transform:translateX(0)} 25%{transform:translateX(-5px)} 75%{transform:translateX(5px)} }
</style></head>
<body>
  <form method="post" autocomplete="off"${wrong ? ' class="err"' : ""}>
    <input name="password" type="password" autofocus required placeholder="Password" aria-label="Password">
    <button type="submit">Enter</button>
  </form>
</body></html>`;
}
