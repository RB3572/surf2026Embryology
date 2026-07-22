// Single source of truth for the site's logins, shared by the edge middleware (the gate) and
// every /api function — so the token→identity map can never drift out of sync.
//
// Two tiers:
//   BASE  — hardcoded accounts. Always work, even if the database is unreachable (fail-safe;
//           an admin can always get in). Resolved WITHOUT any DB call.
//   users — a Neon table the admin can add to at runtime from the console. Cached (30 s) so the
//           gate stays fast; a DB failure just falls back to BASE (never locks anyone out).
//
// Passwords are stored/served in the clear (admin-only) on purpose: this is a small shared-password
// research gate, and the admin needs to read them back to hand them out. Not per-user secret auth.
import { neon } from "@neondatabase/serverless";

export const COOKIE = "surf_gate";
export const MAX_AGE = 60 * 60 * 24 * 30; // remember a login for 30 days

export const BASE = [
  { pw: "Sonichedgehog1", token: "s1_kathytam_b7f3a92c8d1e4056a19d", user: "Kathy Tam", role: "user" },
  { pw: "AdminPass",      token: "adm_owner_5c1e9a2f7b3d8064c2a7f1", user: "Admin",     role: "admin" },
  { pw: "HarryW1!",       token: "hw_harryw_3f9c1e7a5b2d8046e1c9",   user: "Harry W",   role: "user" },
];

const CONN = process.env.DATABASE_URL || process.env.POSTGRES_URL ||
             process.env.DATABASE_URL_UNPOOLED || process.env.POSTGRES_URL_NON_POOLING;
const sql = CONN ? neon(CONN) : null;
export const hasDb = !!sql;

export function cookieVal(header, name) {
  for (const c of (header || "").split(/;\s*/)) {
    const i = c.indexOf("=");
    if (i > 0 && c.slice(0, i) === name) return c.slice(i + 1);
  }
  return null;
}

// ---- dynamically-added users (Neon `users` table) ----
let _cache = null, _at = 0;
export async function dbUsers(force = false) {
  const now = Date.now();
  if (sql && (force || !_cache || now - _at > 30000)) {
    try {
      const rows = await sql`SELECT usr, pw, token, role FROM users ORDER BY created`;
      _cache = rows.map((r) => ({ user: r.usr, pw: r.pw, token: r.token, role: r.role }));
      _at = now;
    } catch (_) { if (!_cache) _cache = []; }   // table missing / DB down → BASE still works
  } else if (!sql && !_cache) _cache = [];
  return _cache || [];
}

export async function accountByToken(token) {
  if (!token) return null;
  const b = BASE.find((a) => a.token === token); if (b) return b;   // BASE resolves with no DB call
  return (await dbUsers()).find((a) => a.token === token) || null;
}
export async function accountByPassword(pw) {
  const b = BASE.find((a) => a.pw === pw); if (b) return b;
  return (await dbUsers(true)).find((a) => a.pw === pw) || null;    // fresh so a just-added login works at once
}
export async function allAccounts() { return BASE.concat(await dbUsers()); }

// ---- admin mutations (from api/users.mjs) ----
async function ensureUsersTable() {
  if (!sql) throw new Error("no database");
  await sql`CREATE TABLE IF NOT EXISTS users (
    usr TEXT PRIMARY KEY, pw TEXT NOT NULL, token TEXT NOT NULL,
    role TEXT NOT NULL DEFAULT 'user', created TIMESTAMPTZ DEFAULT now())`;
}
function genToken(usr) {
  const slug = String(usr || "").toLowerCase().replace(/[^a-z0-9]+/g, "").slice(0, 14) || "user";
  const rnd = (globalThis.crypto && crypto.randomUUID)
    ? crypto.randomUUID().replace(/-/g, "")
    : (Date.now().toString(16) + Math.random().toString(16).slice(2)).padEnd(20, "0");
  return `u_${slug}_${rnd.slice(0, 20)}`;
}
export async function addUser(usr, pw, role) {
  await ensureUsersTable();
  const token = genToken(usr);
  await sql`INSERT INTO users (usr, pw, token, role) VALUES (${usr}, ${pw}, ${token}, ${role})
            ON CONFLICT (usr) DO UPDATE SET pw = ${pw}, role = ${role}`;   // keep the existing token on re-add
  _cache = null;                                                          // invalidate so the login is picked up
  return token;
}
export async function removeUser(usr) {
  await ensureUsersTable();
  await sql`DELETE FROM users WHERE usr = ${usr}`;
  _cache = null;
}
