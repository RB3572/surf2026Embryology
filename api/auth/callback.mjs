// GET /api/auth/callback?code=...&state=...
//
// Where Supabase returns after Google/Apple sign-in. We:
//   1. verify the PKCE/state cookie (CSRF + replay protection)
//   2. exchange the auth code for a Supabase session   (server-side)
//   3. read the user id from that session
//   4. ask PostgREST whether the user is a member of an allowed lab — the query runs AS THE
//      USER, so lab_members' RLS is what actually decides. No service-role key is ever used
//      (that would bypass RLS entirely), and a non-member simply gets zero rows back.
//   5. mint our own httpOnly session cookie (user id + display identity + admin flag)
//
// The Supabase access/refresh tokens are deliberately discarded after step 4: this site is
// read-only and never acts on the user's behalf in Lab Logger.

import {
  PKCE_COOKIE, SESSION_COOKIE, PROVIDER_COOKIE, AUTOTRY_COOKIE,
  verifyToken, signToken, readCookie, cookieHeader, clearCookieHeader,
  requiredEnv, allowedLabIds, safeNext, normalizeProvider, isAdminUser, isMemberResult,
} from "../../lib/session.mjs";
import { bindAccessRow } from "../../lib/projects.mjs";

// How long a session lasts before we silently re-check lab membership. Longer = fewer
// bounces; shorter = removal from the lab takes effect sooner. Re-auth is silent once a
// provider is remembered, so 12h costs the user nothing.
const SESSION_TTL = 60 * 60 * Math.min(720,
  Math.max(1, parseInt(process.env.SESSION_TTL_HOURS, 10) || 12));

const PAGE_CSS = `body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;
background:#f7f8fa;color:#1f2430;display:flex;min-height:100vh;align-items:center;
justify-content:center;margin:0}.c{max-width:28rem;padding:2rem;background:#fff;
border:1px solid #e1e6ed;border-radius:12px;text-align:center}
h1{font-size:1.12rem;margin:0 0 .6rem}p{font-size:.9rem;color:#5b6472;line-height:1.5;margin:0 0 1rem}
a{color:#2f6df3}`;

function fail(res, status, message, detail) {
  res.statusCode = status;
  res.setHeader("content-type", "text/html; charset=utf-8");
  res.setHeader("Cache-Control", "no-store");
  // A silent attempt just failed — forget the remembered provider and the breaker so the
  // user gets the chooser rather than a redirect loop.
  res.setHeader("Set-Cookie", [
    clearCookieHeader(PKCE_COOKIE),
    clearCookieHeader(PROVIDER_COOKIE),
    clearCookieHeader(AUTOTRY_COOKIE),
  ]);
  res.end(`<!doctype html><meta charset="utf-8"><title>Sign-in problem</title>
<style>${PAGE_CSS}</style>
<div class="c"><h1>${message}</h1>${detail ? `<p>${escapeHtml(detail)}</p>` : ""}
<p><a href="/login">Try again</a></p></div>`);
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, (c) => (
    { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
}

export default async function handler(req, res) {
  const missing = requiredEnv();
  if (missing.length) return fail(res, 500, "Auth is not configured", missing.join(", "));

  const secret = process.env.SESSION_SECRET;
  const SUPABASE_URL = process.env.SUPABASE_URL;
  const ANON = process.env.SUPABASE_ANON_KEY;

  // 1. PKCE + CSRF state
  const pkce = await verifyToken(readCookie(req.headers.cookie, PKCE_COOKIE), secret);
  if (!pkce) return fail(res, 400, "Your sign-in link expired", "Please start again.");
  if (!req.query.state || req.query.state !== pkce.s) {
    return fail(res, 400, "Sign-in could not be verified", "State mismatch.");
  }
  if (req.query.error) {
    return fail(res, 400, "Sign-in was cancelled",
      String(req.query.error_description || req.query.error));
  }
  const code = String(req.query.code || "");
  if (!code) return fail(res, 400, "Sign-in failed", "No authorization code returned.");

  // 2. Exchange the Supabase PKCE code for a session.
  let session;
  try {
    const r = await fetch(`${SUPABASE_URL}/auth/v1/token?grant_type=pkce`, {
      method: "POST",
      headers: { "Content-Type": "application/json", apikey: ANON },
      body: JSON.stringify({ auth_code: code, code_verifier: pkce.v }),
    });
    if (!r.ok) {
      return fail(res, 502, "Sign-in failed", `Supabase rejected the exchange (${r.status}).`);
    }
    session = await r.json();
  } catch (e) {
    return fail(res, 502, "Sign-in failed", e?.message || String(e));
  }
  const accessToken = session?.access_token;
  if (!accessToken) return fail(res, 502, "Sign-in failed", "No session returned.");

  // 3. Who is this? Name/email come straight from the provider profile.
  const userId = session?.user?.id;
  if (!userId) return fail(res, 502, "Sign-in failed", "No user id in session.");
  const meta = session.user.user_metadata || {};
  const email = session.user.email || "";
  const name = meta.full_name || meta.name || (email ? email.split("@")[0] : "member");

  // 4. Membership check, executed under the user's own RLS.
  const labs = allowedLabIds();
  let isMember = false;
  let labRole = "";
  try {
    const q = new URL(`${SUPABASE_URL}/rest/v1/lab_members`);
    q.searchParams.set("select", "lab_id,role");
    q.searchParams.set("lab_id", `in.(${labs.join(",")})`);
    q.searchParams.set("limit", "1");
    const r = await fetch(q, {
      headers: { apikey: ANON, Authorization: `Bearer ${accessToken}` },
    });
    if (!r.ok) {
      return fail(res, 502, "Could not verify lab access", `Lookup failed (${r.status}).`);
    }
    const rows = await r.json();
    // RLS returns rows ONLY for labs this user actually belongs to, so a non-empty result is
    // the membership proof. An empty array means NOT a member — see isMemberResult.
    isMember = isMemberResult(rows);
    if (isMember) labRole = String(rows[0].role || "");
  } catch (e) {
    return fail(res, 502, "Could not verify lab access", e?.message || String(e));
  }

  if (!isMember) {
    res.statusCode = 403;
    res.setHeader("content-type", "text/html; charset=utf-8");
    res.setHeader("Cache-Control", "no-store");
    // No session cookie is minted. Clear the silent-sign-in helpers too.
    res.setHeader("Set-Cookie", [
      clearCookieHeader(PKCE_COOKIE),
      clearCookieHeader(PROVIDER_COOKIE),
      clearCookieHeader(AUTOTRY_COOKIE),
    ]);
    res.end(`<!doctype html><meta charset="utf-8"><title>No access</title>
<style>${PAGE_CSS}</style>
<div class="c"><h1>You don't have access to this site</h1>
<p>This site is limited to members of the Zernicka-Goetz Lab in
<a href="https://www.lab-logger.com">Lab&nbsp;Logger</a>. Ask your PI or a lab admin to add
you to the lab, then sign in again.</p>
<p>Signed in as ${escapeHtml(email || name)}.</p>
<p><a href="/login?signedout=1">Back to sign in</a></p></div>`);
    return;
  }

  // 5. Our own session cookie. Carries the display identity (so the admin console is readable
  //    without re-querying Lab Logger) and an admin flag. Admin is the ADMIN_USER_IDS
  //    allowlist ONLY — a Lab Logger admin/PI role does not confer it here.
  const isAdmin = isAdminUser(userId);

  // Bind any email-seeded project-access row to this user id, so a restriction created before
  // the person ever signed in applies from their very first visit. Never fatal.
  try { await bindAccessRow(userId, email, name); } catch (_) { /* access is best-effort */ }

  const token = await signToken({
    sub: userId, em: email, nm: name, rl: labRole, adm: isAdmin,
    exp: Math.floor(Date.now() / 1000) + SESSION_TTL,
  }, secret);

  const usedProvider = normalizeProvider(pkce.p);
  res.setHeader("Set-Cookie", [
    cookieHeader(SESSION_COOKIE, token, SESSION_TTL),
    clearCookieHeader(PKCE_COOKIE),
    // Remember the provider so future visits skip the chooser entirely.
    ...(usedProvider ? [cookieHeader(PROVIDER_COOKIE, usedProvider, 60 * 60 * 24 * 365)] : []),
    clearCookieHeader(AUTOTRY_COOKIE),
  ]);
  res.setHeader("Cache-Control", "no-store");
  res.statusCode = 302;
  res.setHeader("Location", safeNext(pkce.n));
  res.end();
}
