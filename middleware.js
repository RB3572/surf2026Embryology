// Password-gate the whole site at the edge (Vercel Routing/Edge Middleware).
//
// This runs BEFORE any file is served, so nothing — pages, scripts, styles, or the
// data files under /data — is reachable without the shared password. It is framework-
// agnostic: plain Web APIs, no imports. Returning nothing lets an unlocked request
// continue to the static file; returning a Response short-circuits (the gate / redirect).
//
// The password lives only in this edge module (server-side); it is never shipped to the
// browser. To change or remove protection, edit PASSWORD below (or delete this file).

const PASSWORD = "SonicHedgehog";
const COOKIE = "surf_gate";
const TOKEN = "unlocked-7b3f9a2c8d1e4056"; // opaque unlock marker; stays on the edge
const MAX_AGE = 60 * 60 * 24 * 30;         // remember an unlocked visitor for 30 days

export const config = {
  // gate every path except Vercel's own internals
  matcher: "/((?!_vercel).*)",
};

export default async function middleware(request) {
  const url = new URL(request.url);
  const cookies = request.headers.get("cookie") || "";
  const unlocked = cookies.split(/;\s*/).includes(`${COOKIE}=${TOKEN}`);

  if (unlocked) return; // ✓ already unlocked → continue, serve the requested file

  // login form submitted
  if (request.method === "POST") {
    let pwd = "";
    try {
      const form = await request.formData();
      pwd = String(form.get("password") || "");
    } catch (_) {
      /* ignore malformed body */
    }
    if (pwd === PASSWORD) {
      return new Response(null, {
        status: 303,
        headers: {
          "set-cookie": `${COOKIE}=${TOKEN}; Path=/; Max-Age=${MAX_AGE}; HttpOnly; Secure; SameSite=Lax`,
          location: url.pathname + url.search,
          "cache-control": "no-store",
        },
      });
    }
    return gate(true); // wrong password
  }

  return gate(false); // not unlocked → present the gate immediately
}

function gate(wrong) {
  return new Response(gateHtml(wrong), {
    status: 401,
    headers: {
      "content-type": "text/html; charset=utf-8",
      "cache-control": "no-store",
    },
  });
}

function gateHtml(wrong) {
  const err = wrong
    ? `<p class="err" role="alert">Incorrect password — try again.</p>`
    : "";
  return `<!doctype html><html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Private · SURF 2026 · Spatial Transcriptomics</title>
<style>
  :root { color-scheme: dark; }
  * { box-sizing: border-box; }
  html, body { height: 100%; margin: 0; }
  body {
    font: 15px/1.5 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    color: #e7e9ee; display: grid; place-items: center; padding: 24px;
    background: radial-gradient(1200px 600px at 50% -10%, #1b2440 0%, #0b1020 55%, #070a14 100%);
  }
  .card {
    width: 100%; max-width: 380px; background: rgba(21,27,46,.72);
    border: 1px solid rgba(255,255,255,.08); border-radius: 16px; padding: 30px 28px 26px;
    box-shadow: 0 20px 60px rgba(0,0,0,.45); backdrop-filter: blur(8px);
  }
  .lock { width: 44px; height: 44px; border-radius: 12px; display: grid; place-items: center;
    background: linear-gradient(140deg, #f97316, #db2777); margin-bottom: 16px; font-size: 22px; }
  h1 { font-size: 18px; margin: 0 0 4px; letter-spacing: .2px; }
  .sub { margin: 0 0 20px; color: #9aa3b8; font-size: 13px; }
  label { display: block; font-size: 12px; color: #9aa3b8; margin: 0 0 6px; }
  input[type=password] {
    width: 100%; padding: 11px 13px; border-radius: 10px; font-size: 15px;
    background: #0d1324; border: 1px solid rgba(255,255,255,.12); color: #fff; outline: none;
  }
  input[type=password]:focus { border-color: #f97316; box-shadow: 0 0 0 3px rgba(249,115,22,.22); }
  button {
    width: 100%; margin-top: 14px; padding: 11px 14px; border: 0; border-radius: 10px; cursor: pointer;
    font-size: 15px; font-weight: 600; color: #fff; background: linear-gradient(140deg, #f97316, #db2777);
  }
  button:hover { filter: brightness(1.06); }
  .err { color: #fca5a5; font-size: 13px; margin: 12px 0 0; }
  .foot { margin: 18px 0 0; color: #5f677d; font-size: 11px; text-align: center; }
</style></head>
<body>
  <form class="card" method="post" autocomplete="off">
    <div class="lock" aria-hidden="true">🔒</div>
    <h1>This site is private</h1>
    <p class="sub">SURF 2026 · Spatial Transcriptomics. Enter the password to continue.</p>
    <label for="pw">Password</label>
    <input id="pw" name="password" type="password" autofocus required aria-label="Password">
    <button type="submit">Unlock</button>
    ${err}
    <p class="foot">Access is restricted to the research team.</p>
  </form>
</body></html>`;
}
