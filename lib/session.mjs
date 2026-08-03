// Shared auth helpers for the SURF 2026 embryology site.
//
// Access is delegated to Lab Logger's Supabase project: you are allowed in if you are a
// member of an allowed lab, and nothing else. There is no user list and no password here.
//
// Everything uses Web Crypto (globalThis.crypto.subtle) rather than node:crypto, so the SAME
// module works in both runtimes we need it in: the Vercel Node functions under api/ and the
// Edge middleware at the root.
//
// The session cookie is a compact signed token:  base64url(JSON).base64url(HMAC)
// It carries the Supabase user id + expiry, plus display name/email/lab role for the admin
// console. It never carries a Supabase access or refresh token.

const enc = new TextEncoder();
const dec = new TextDecoder();

export const SESSION_COOKIE = "surf_session";
export const PKCE_COOKIE = "surf_pkce";
/** Which provider this browser last signed in with, so we can send them straight there next
 *  time instead of showing a chooser. Not sensitive. */
export const PROVIDER_COOKIE = "surf_provider";
/** Set for ~60s while a silent sign-in is in flight; if we come back round still
 *  unauthenticated it stops us from looping forever. */
export const AUTOTRY_COOKIE = "surf_auto";

// ── base64url ──────────────────────────────────────────────────────────────
function b64urlFromBytes(input) {
  const bytes = input instanceof Uint8Array ? input : new Uint8Array(input);
  let s = "";
  for (const b of bytes) s += String.fromCharCode(b);
  return btoa(s).replace(/\+/g, "-").replace(/\//g, "_").replace(/=+$/, "");
}

function bytesFromB64url(s) {
  const pad = s.length % 4 ? "=".repeat(4 - (s.length % 4)) : "";
  const bin = atob(s.replace(/-/g, "+").replace(/_/g, "/") + pad);
  const out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out;
}

// ── HMAC-signed tokens ─────────────────────────────────────────────────────
async function hmacKey(secret) {
  return crypto.subtle.importKey(
    "raw", enc.encode(secret),
    { name: "HMAC", hash: "SHA-256" }, false, ["sign", "verify"]);
}

/** Sign a payload object. `exp` (unix seconds) is required by verifyToken(). */
export async function signToken(payload, secret) {
  const body = b64urlFromBytes(enc.encode(JSON.stringify(payload)));
  const sig = await crypto.subtle.sign("HMAC", await hmacKey(secret), enc.encode(body));
  return `${body}.${b64urlFromBytes(sig)}`;
}

/** Verify + decode. Returns the payload, or null if tampered/expired/malformed. */
export async function verifyToken(token, secret) {
  if (!token || typeof token !== "string" || !secret) return null;
  const dot = token.indexOf(".");
  if (dot < 1) return null;
  const body = token.slice(0, dot);
  try {
    const ok = await crypto.subtle.verify(
      "HMAC", await hmacKey(secret),
      bytesFromB64url(token.slice(dot + 1)), enc.encode(body));
    if (!ok) return null;                       // constant-time inside subtle
    const payload = JSON.parse(dec.decode(bytesFromB64url(body)));
    if (typeof payload?.exp !== "number") return null;
    if (payload.exp < Math.floor(Date.now() / 1000)) return null;   // expired
    return payload;
  } catch {
    return null;
  }
}

// ── PKCE ───────────────────────────────────────────────────────────────────
export function randomToken(bytes = 32) {
  return b64urlFromBytes(crypto.getRandomValues(new Uint8Array(bytes)));
}

export async function sha256b64url(text) {
  return b64urlFromBytes(await crypto.subtle.digest("SHA-256", enc.encode(text)));
}

// ── cookies ────────────────────────────────────────────────────────────────
export function readCookie(header, name) {
  if (!header) return null;
  for (const part of header.split(";")) {
    const i = part.indexOf("=");
    if (i < 0) continue;
    if (part.slice(0, i).trim() === name) {
      return decodeURIComponent(part.slice(i + 1).trim());
    }
  }
  return null;
}

/** httpOnly + Secure + SameSite=Lax. Lax (not Strict) so the cookie is sent on the top-level
 *  redirect back from Supabase; it still blocks cross-site POSTs. */
export function cookieHeader(name, value, maxAgeSec) {
  return [
    `${name}=${encodeURIComponent(value)}`,
    "Path=/", "HttpOnly", "Secure", "SameSite=Lax",
    `Max-Age=${maxAgeSec}`,
  ].join("; ");
}

export function clearCookieHeader(name) {
  return `${name}=; Path=/; HttpOnly; Secure; SameSite=Lax; Max-Age=0`;
}

/** google | apple, or null. Used to decide whether we can skip the chooser.
 *  NEVER guess: sending an Apple user to Google can silently authenticate a DIFFERENT
 *  account that isn't in the lab. No signal => return null => show the chooser. */
export function normalizeProvider(v) {
  return v === "google" || v === "apple" ? v : null;
}

// ── env ────────────────────────────────────────────────────────────────────
export function requiredEnv() {
  const missing = [];
  for (const k of ["SUPABASE_URL", "SUPABASE_ANON_KEY", "SESSION_SECRET", "ALLOWED_LAB_IDS"]) {
    if (!process.env[k]) missing.push(k);
  }
  return missing;
}

export function allowedLabIds() {
  return (process.env.ALLOWED_LAB_IDS || "")
    .split(",").map((s) => s.trim()).filter(Boolean);
}

/** Public origin of this deployment. */
export function originFor(req) {
  if (process.env.SITE_ORIGIN) return process.env.SITE_ORIGIN;
  const host = req.headers["x-forwarded-host"] || req.headers.host || "incrementum.rishib.com";
  const proto = req.headers["x-forwarded-proto"] || "https";
  return `${proto}://${host}`;
}

/** Only same-site relative paths survive as a post-login destination, so ?next= can't be
 *  turned into an open redirect. Rejects "//evil.com" (protocol-relative) and "https://…". */
export function safeNext(next) {
  if (typeof next !== "string") return "/";
  if (!next.startsWith("/") || next.startsWith("//")) return "/";
  return next;
}

/** Is this user an admin?
 *
 *  ONLY the explicit ADMIN_USER_IDS allowlist counts. A Lab Logger lab role of admin/PI does
 *  NOT confer admin here: lab membership is managed by other people, so inheriting it would
 *  silently hand this site's analytics console to anyone promoted in Lab Logger.
 *
 *  Fails closed: with ADMIN_USER_IDS unset, nobody is an admin and /admin 404s for everyone. */
export function isAdminUser(userId) {
  if (!userId) return false;
  return (process.env.ADMIN_USER_IDS || "")
    .split(",").map((s) => s.trim()).filter(Boolean)
    .includes(userId);
}

/** The membership proof. PostgREST runs the lab_members query under the USER'S OWN token, so
 *  RLS returns rows only for labs they actually belong to.
 *
 *  ⚠️ A non-member gets HTTP 200 with an EMPTY ARRAY — not an error. Checking `Array.isArray(rows)`
 *  alone would admit EVERYONE. The length check is the entire access control decision. */
export function isMemberResult(rows) {
  return Array.isArray(rows) && rows.length > 0;
}
