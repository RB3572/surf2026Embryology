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

const COOKIE = "surf_gate";
const MAX_AGE = 60 * 60 * 24 * 30; // remember a login for 30 days

const ACCOUNTS = [
  { pw: "Sonichedgehog1", token: "s1_kathytam_b7f3a92c8d1e4056a19d", user: "Kathy Tam", role: "user" },
  { pw: "AdminPass",      token: "adm_owner_5c1e9a2f7b3d8064c2a7f1", user: "Admin",     role: "admin" },
];
const BY_TOKEN = new Map(ACCOUNTS.map((a) => [a.token, a]));

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

  // Admin-only surface. Anything under /admin (page, assets, data) is 404 unless this is an
  // admin session — indistinguishable from "does not exist", so non-admins never see it.
  const p = url.pathname;
  const adminOnly = p === "/admin" || p.startsWith("/admin/") || p.startsWith("/admin.");
  if (adminOnly && (!account || account.role !== "admin")) {
    return new Response("Not Found", {
      status: 404,
      headers: { "content-type": "text/plain; charset=utf-8", "cache-control": "no-store" },
    });
  }

  if (account) return; // authenticated → serve the requested file

  // login form submitted
  if (request.method === "POST") {
    let pwd = "";
    try { pwd = String((await request.formData()).get("password") || ""); } catch (_) { /* malformed */ }
    const match = ACCOUNTS.find((a) => a.pw === pwd);
    if (match) {
      return new Response(null, {
        status: 303,
        headers: {
          "set-cookie": `${COOKIE}=${match.token}; Path=/; Max-Age=${MAX_AGE}; HttpOnly; Secure; SameSite=Lax`,
          location: url.pathname + url.search,
          "cache-control": "no-store",
        },
      });
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
