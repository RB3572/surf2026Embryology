// Enforcement tests for the edge middleware — the component that actually decides who sees
// the dataset. These drive the real exported middleware function with real Request objects.
//
// DATABASE_URL is deliberately left unset here, so the per-project access lookup short-circuits
// to "no restriction" and these tests exercise the authentication path only.

import test from "node:test";
import assert from "node:assert/strict";

process.env.SESSION_SECRET = "test-secret-do-not-use-in-production";
const SECRET = process.env.SESSION_SECRET;

const { signToken } = await import("../lib/session.mjs");
const { default: middleware } = await import("../middleware.js");

const future = () => Math.floor(Date.now() / 1000) + 3600;

const req = (path, { accept = "", cookie = "" } = {}) =>
  new Request(`https://incrementum.rishib.com${path}`, {
    headers: { ...(accept ? { accept } : {}), ...(cookie ? { cookie } : {}) },
  });

const HTML = "text/html,application/xhtml+xml";

async function sessionCookie(payload = {}) {
  const t = await signToken({ sub: "user-1", nm: "Test", adm: false, exp: future(), ...payload }, SECRET);
  return `surf_session=${encodeURIComponent(t)}`;
}

// ── fail closed ────────────────────────────────────────────────────────────
test("returns 503 when SESSION_SECRET is unset, rather than serving data", async () => {
  const saved = process.env.SESSION_SECRET;
  delete process.env.SESSION_SECRET;
  try {
    const res = await middleware(req("/data/zygote_sperm.json"));
    assert.equal(res.status, 503);
  } finally { process.env.SESSION_SECRET = saved; }
});

// ── unauthenticated ────────────────────────────────────────────────────────
test("a data file is 401 for an unauthenticated fetch — not 200", async () => {
  // This is the single most important assertion in the suite: curl must not get the dataset.
  const res = await middleware(req("/data/zygote_sperm.json"));
  assert.equal(res.status, 401);
});

test("assets and project pages are 401 without an HTML Accept header", async () => {
  for (const p of ["/viewer-core.js", "/landing.css", "/sperm-sphere.html", "/data/x.json.gz"]) {
    const res = await middleware(req(p));
    assert.equal(res.status, 401, `${p} should be 401`);
  }
});

test("an HTML navigation redirects to /login carrying the original path", async () => {
  const res = await middleware(req("/sperm-sphere.html?gene=Actb", { accept: HTML }));
  assert.equal(res.status, 302);
  const loc = new URL(res.headers.get("location"));
  assert.equal(loc.pathname, "/login");
  assert.equal(loc.searchParams.get("next"), "/sperm-sphere.html?gene=Actb");
});

test("the admin page does not leak its existence to anonymous callers", async () => {
  const res = await middleware(req("/admin", { accept: HTML }));
  assert.equal(res.status, 302);           // same treatment as any other page
  assert.match(res.headers.get("location"), /\/login/);
});

// ── silent continuation ────────────────────────────────────────────────────
test("a remembered provider continues silently instead of showing the chooser", async () => {
  const res = await middleware(req("/", { accept: HTML, cookie: "surf_provider=google" }));
  assert.equal(res.status, 302);
  const loc = new URL(res.headers.get("location"));
  assert.equal(loc.pathname, "/api/auth/login");
  assert.equal(loc.searchParams.get("provider"), "google");
  assert.match(res.headers.get("set-cookie") || "", /surf_auto=1/);   // breaker armed
});

test("a ?provider= hint on the incoming link is honoured", async () => {
  const res = await middleware(req("/?provider=apple", { accept: HTML }));
  const loc = new URL(res.headers.get("location"));
  assert.equal(loc.pathname, "/api/auth/login");
  assert.equal(loc.searchParams.get("provider"), "apple");
});

test("a bogus ?provider= is never guessed — it falls back to the chooser", async () => {
  const res = await middleware(req("/?provider=facebook", { accept: HTML }));
  assert.match(new URL(res.headers.get("location")).pathname, /^\/login$/);
});

test("the auto-try breaker stops a redirect loop after a failed silent attempt", async () => {
  const res = await middleware(req("/", {
    accept: HTML, cookie: "surf_provider=google; surf_auto=1",
  }));
  assert.equal(new URL(res.headers.get("location")).pathname, "/login");
});

test("signing out is possible — ?signedout=1 is never auto-continued", async () => {
  // Without this guard the remembered provider would bounce the user straight back in.
  const res = await middleware(req("/login?signedout=1", {
    accept: HTML, cookie: "surf_provider=google",
  }));
  assert.notEqual(res.status, 302);   // serves the login page instead of re-authenticating
});

test("/login is served, never redirected to itself", async () => {
  const res = await middleware(req("/login", { accept: HTML }));
  assert.notEqual(res.status, 302);
  assert.equal(res.headers.get("cache-control"), "no-store");
});

test("/login?next= cannot be used as an open redirect", async () => {
  const res = await middleware(req("/login?next=//evil.com", {
    accept: HTML, cookie: "surf_provider=google",
  }));
  const loc = new URL(res.headers.get("location"));
  assert.equal(loc.searchParams.get("next"), "/");   // sanitised, not //evil.com
});

// ── authenticated ──────────────────────────────────────────────────────────
test("a valid session passes through and is marked Vary: Cookie", async () => {
  const res = await middleware(req("/data/zygote_sperm.json", { cookie: await sessionCookie() }));
  assert.notEqual(res.status, 401);
  assert.equal(res.headers.get("vary"), "Cookie");
});

test("an expired session is treated as anonymous", async () => {
  const stale = await signToken({ sub: "u", exp: Math.floor(Date.now() / 1000) - 10 }, SECRET);
  const res = await middleware(req("/data/x.json", { cookie: `surf_session=${stale}` }));
  assert.equal(res.status, 401);
});

test("a session signed with the wrong secret is treated as anonymous", async () => {
  const forged = await signToken({ sub: "u", adm: true, exp: future() }, "attacker-secret");
  const res = await middleware(req("/data/x.json", { cookie: `surf_session=${forged}` }));
  assert.equal(res.status, 401);
});

// ── admin surface ──────────────────────────────────────────────────────────
test("a non-admin member gets 404 (not 403) on the admin page and API", async () => {
  const cookie = await sessionCookie({ adm: false });
  for (const p of ["/admin", "/admin.html", "/admin.card.html", "/api/admin"]) {
    const res = await middleware(req(p, { accept: HTML, cookie }));
    assert.equal(res.status, 404, `${p} should 404 for a non-admin`);
  }
});

test("an admin reaches the admin surface", async () => {
  const cookie = await sessionCookie({ adm: true });
  for (const p of ["/admin", "/admin.html", "/api/admin"]) {
    const res = await middleware(req(p, { accept: HTML, cookie }));
    assert.notEqual(res.status, 404, `${p} should be reachable by an admin`);
  }
});
