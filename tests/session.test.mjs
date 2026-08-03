// Unit tests for the auth primitives. These cover the three things that, if wrong, silently
// let the wrong people in:
//   1. session token forgery / expiry
//   2. the ?next= open-redirect guard
//   3. the lab-membership check — where an empty array MUST mean "denied"
//
// Run with:  npm test        (node --test tests/)

import test from "node:test";
import assert from "node:assert/strict";

import {
  signToken, verifyToken, safeNext, normalizeProvider,
  readCookie, cookieHeader, clearCookieHeader, isMemberResult, isAdminUser,
} from "../lib/session.mjs";

const SECRET = "test-secret-do-not-use-in-production";
const future = () => Math.floor(Date.now() / 1000) + 3600;
const past = () => Math.floor(Date.now() / 1000) - 60;

// ── session sign / verify ──────────────────────────────────────────────────
test("round-trips a payload", async () => {
  const t = await signToken({ sub: "user-1", adm: false, exp: future() }, SECRET);
  const p = await verifyToken(t, SECRET);
  assert.equal(p.sub, "user-1");
  assert.equal(p.adm, false);
});

test("rejects a tampered payload", async () => {
  const t = await signToken({ sub: "user-1", adm: false, exp: future() }, SECRET);
  // Forge an admin session by swapping the body but keeping the signature.
  const forgedBody = Buffer.from(JSON.stringify({ sub: "user-1", adm: true, exp: future() }))
    .toString("base64url");
  const forged = `${forgedBody}.${t.slice(t.indexOf(".") + 1)}`;
  assert.equal(await verifyToken(forged, SECRET), null);
});

test("rejects a tampered signature", async () => {
  const t = await signToken({ sub: "user-1", exp: future() }, SECRET);
  const [body, sig] = t.split(".");
  const flipped = sig.slice(0, -1) + (sig.endsWith("A") ? "B" : "A");
  assert.equal(await verifyToken(`${body}.${flipped}`, SECRET), null);
});

test("rejects a token signed with a different secret", async () => {
  const t = await signToken({ sub: "user-1", exp: future() }, "some-other-secret");
  assert.equal(await verifyToken(t, SECRET), null);
});

test("rejects an expired token", async () => {
  const t = await signToken({ sub: "user-1", exp: past() }, SECRET);
  assert.equal(await verifyToken(t, SECRET), null);
});

test("rejects a token with no exp — a session must always expire", async () => {
  const t = await signToken({ sub: "user-1" }, SECRET);
  assert.equal(await verifyToken(t, SECRET), null);
});

test("rejects a non-numeric exp", async () => {
  const t = await signToken({ sub: "user-1", exp: "9999999999" }, SECRET);
  assert.equal(await verifyToken(t, SECRET), null);
});

test("rejects malformed / missing / empty tokens", async () => {
  for (const bad of [null, undefined, "", "no-dot", ".", "a.b.c", 42, {}]) {
    assert.equal(await verifyToken(bad, SECRET), null, `should reject ${JSON.stringify(bad)}`);
  }
});

test("rejects any token when the secret is missing (fail closed)", async () => {
  const t = await signToken({ sub: "user-1", exp: future() }, SECRET);
  assert.equal(await verifyToken(t, ""), null);
  assert.equal(await verifyToken(t, undefined), null);
});

// ── ?next= open-redirect guard ─────────────────────────────────────────────
test("safeNext keeps same-site paths", () => {
  assert.equal(safeNext("/"), "/");
  assert.equal(safeNext("/sperm-sphere.html"), "/sperm-sphere.html");
  assert.equal(safeNext("/page?x=y&z=1"), "/page?x=y&z=1");
  assert.equal(safeNext("/a/b/c#frag"), "/a/b/c#frag");
});

test("safeNext rejects off-site destinations", () => {
  for (const bad of [
    "//evil.com",                 // protocol-relative — the classic bypass
    "///evil.com",
    "https://evil.com",
    "http://evil.com",
    "javascript:alert(1)",
    "evil.com",
    "\\\\evil.com",
    "",
  ]) {
    assert.equal(safeNext(bad), "/", `should reject ${JSON.stringify(bad)}`);
  }
});

test("safeNext rejects non-strings", () => {
  for (const bad of [null, undefined, 42, {}, ["/ok"]]) assert.equal(safeNext(bad), "/");
});

// ── membership check — the branch that matters most ────────────────────────
test("isMemberResult DENIES an empty array", () => {
  // PostgREST answers 200 + [] for a non-member. Checking Array.isArray(rows) alone here
  // would admit EVERYONE — this assertion is the guard against that regression.
  assert.equal(isMemberResult([]), false);
});

test("isMemberResult ALLOWS a non-empty array", () => {
  assert.equal(isMemberResult([{ lab_id: "3285888c-7cd3-44c4-bb94-7b544d2a645b", role: "member" }]), true);
});

test("isMemberResult denies anything that isn't an array", () => {
  for (const bad of [null, undefined, {}, "rows", 1, true, { rows: [{}] },
                     { message: "permission denied" }]) {
    assert.equal(isMemberResult(bad), false, `should deny ${JSON.stringify(bad)}`);
  }
});

// ── provider handling — never guess ────────────────────────────────────────
test("normalizeProvider accepts only the two real providers", () => {
  assert.equal(normalizeProvider("google"), "google");
  assert.equal(normalizeProvider("apple"), "apple");
});

test("normalizeProvider returns null for anything else, so we show the chooser", () => {
  // Guessing a provider can silently authenticate a DIFFERENT account that isn't in the lab.
  for (const bad of ["Google", "GOOGLE", "facebook", "", null, undefined, 1, {}]) {
    assert.equal(normalizeProvider(bad), null, `should not infer from ${JSON.stringify(bad)}`);
  }
});

// ── admin derivation ───────────────────────────────────────────────────────
test("only the ADMIN_USER_IDS allowlist confers admin", () => {
  const prev = process.env.ADMIN_USER_IDS;
  process.env.ADMIN_USER_IDS = "u-alpha, u-beta";
  assert.equal(isAdminUser("u-alpha"), true);
  assert.equal(isAdminUser("u-beta"), true);      // whitespace around the comma is tolerated
  assert.equal(isAdminUser("u-gamma"), false);
  process.env.ADMIN_USER_IDS = prev;
});

test("a Lab Logger admin/PI role does NOT confer admin here", () => {
  // Lab membership is managed by other people. If a lab role granted this console, promoting
  // somebody in Lab Logger would silently hand them the site's analytics.
  const prev = process.env.ADMIN_USER_IDS;
  process.env.ADMIN_USER_IDS = "the-only-admin";
  for (const role of ["admin", "pi", "owner", "member", ""]) {
    assert.equal(isAdminUser("someone-else", role), false,
      `lab role "${role}" must not confer admin`);
  }
  assert.equal(isAdminUser("the-only-admin", "member"), true);
  process.env.ADMIN_USER_IDS = prev;
});

test("nobody is admin when ADMIN_USER_IDS is unset (fail closed)", () => {
  const prev = process.env.ADMIN_USER_IDS;
  delete process.env.ADMIN_USER_IDS;
  assert.equal(isAdminUser("anyone"), false);
  assert.equal(isAdminUser("anyone", "admin"), false);
  process.env.ADMIN_USER_IDS = prev;
});

test("isAdminUser rejects an empty/missing user id", () => {
  const prev = process.env.ADMIN_USER_IDS;
  // An empty allowlist entry must never match an empty subject.
  process.env.ADMIN_USER_IDS = "a,,b";
  assert.equal(isAdminUser(""), false);
  assert.equal(isAdminUser(undefined), false);
  assert.equal(isAdminUser(null), false);
  process.env.ADMIN_USER_IDS = prev;
});

// ── cookies ────────────────────────────────────────────────────────────────
test("readCookie picks the right value", () => {
  const h = "a=1; surf_session=abc.def; surf_provider=google";
  assert.equal(readCookie(h, "surf_session"), "abc.def");
  assert.equal(readCookie(h, "surf_provider"), "google");
  assert.equal(readCookie(h, "nope"), null);
  assert.equal(readCookie(null, "surf_session"), null);
});

test("readCookie is not fooled by a name that is a suffix of another", () => {
  assert.equal(readCookie("xsurf_session=evil; surf_session=good", "surf_session"), "good");
});

test("session cookies are httpOnly, Secure and SameSite=Lax", () => {
  const c = cookieHeader("surf_session", "v", 60);
  assert.match(c, /HttpOnly/);
  assert.match(c, /Secure/);
  assert.match(c, /SameSite=Lax/);
  assert.match(c, /Max-Age=60/);
  assert.match(c, /Path=\//);
});

test("clearing a cookie expires it immediately", () => {
  assert.match(clearCookieHeader("surf_session"), /Max-Age=0/);
});

test("cookie values are encoded so they cannot inject attributes", () => {
  const c = cookieHeader("surf_session", "a; Domain=evil.com", 60);
  assert.ok(!c.includes("Domain=evil.com"), "value must not break out of the cookie");
});
