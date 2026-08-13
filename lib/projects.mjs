// Per-person project access — the one piece of authorization that is NOT delegated to Lab
// Logger. Lab membership decides whether you get into the site at all; this decides which
// projects a given member may open.
//
// Keyed on the SUPABASE USER ID (not a display name, not an email — Apple's Hide My Email
// means an email can change identity between providers). A member with NO row sees EVERY
// project; a row restricts them to exactly the listed keys.
//
// Pre-seeding: an admin can restrict someone BEFORE they have ever signed in, by email. That
// row is stored under the placeholder key `pending:<lowercased email>` and is rebound to the
// real user id by bindAccessRow() on their first successful sign-in — so the restriction
// applies from their very first visit rather than after a gap.
//
// NOTE: the previous password-era table (`access`, keyed on display name) is intentionally
// left untouched as a historical record. This module uses a new table, `project_access`.

import { neon } from "@neondatabase/serverless";

/** Canonical project list: key = the page basename without .html. Kept in step with the
 *  actual pages on disk; PROJECT_OF below is derived from it so the two cannot drift. */
export const PROJECTS = [
  { key: "pronuclei", label: "Interpronuclei Distance Clock" },
  { key: "extpt", label: "Extended-Time Pronuclei-COM Clock" },
  { key: "pronuclei-assignments", label: "Maternal / Paternal Pronucleus ID" },
  { key: "segments", label: "Gene Enrichment by Segment" },
  { key: "stage-expression", label: "Stage Expression Explorer" },
  { key: "diffusion", label: "mRNA Diffusion Rates" },
  { key: "zygote-planes", label: "Division Plane Sweep (18 candidates)" },
  { key: "planes-all", label: "Division Plane Sweep — every plane" },
  { key: "equatorial-planes", label: "Equatorial Division Plane" },
  { key: "compare-planes", label: "Compare Division Planes" },
  { key: "sperm-division", label: "Sperm-Defined Division Plane" },
  { key: "pseudosperm", label: "Pseudosperm Division Plane" },
  { key: "sperm-sphere", label: "Sperm-Entry-Site Enrichment" },
  { key: "sperm-pseudotime", label: "Sperm Location vs Pseudotime" },
  { key: "sperm-pca", label: "Sperm Location Prediction (PCA)" },
  { key: "alphabeta", label: "2-Cell Blastomere α/β Labelling" },
  { key: "clustering", label: "Spatial Gene Clustering" },
  { key: "contact", label: "Blastomere Contact Enrichment" },
  { key: "sperm-map", label: "Sperm Location Browser" },
  { key: "alignment", label: "Sperm Alignment (2-cell)" },
  { key: "size", label: "Embryo Size" },
  { key: "scheffler", label: "Scheffler Pseudotime Model Training" },
];

export const PROJECT_KEYS = new Set(PROJECTS.map((p) => p.key));

/** Pages that are admin-only outright — not part of the per-user allowlist system at all, the
 *  same tier as /admin itself (middleware 404s them for anyone without session.adm). Kept OUT
 *  of PROJECTS so the admin console's per-user access matrix doesn't offer a pointless toggle
 *  for a page nobody but an admin can ever open. Their cards live in admin.card.html instead of
 *  the public landing groups. */
export const ADMIN_ONLY_PAGES = new Set([
  "pronuclear-pseudotime", "pseudotime-calibration", "vision-pseudotime", "pn3d-transcripts",
  "axes",
  // ---- ported from HighResSlideshowExports, admin-only until reviewed ----
  "clocktx",
]);

/** "<key>.html" -> key, for the middleware's path check. Derived, never hand-maintained. */
export const PROJECT_OF = Object.fromEntries(PROJECTS.map((p) => [`${p.key}.html`, p.key]));

/** Which project (if any) a request path corresponds to. Handles both "/x.html" and the
 *  extensionless "/x" form, so adding cleanUrls later can't silently open a hole. */
export function projectForPath(pathname) {
  const bare = String(pathname || "").replace(/^\/+/, "").replace(/\/+$/, "");
  if (PROJECT_OF[bare]) return PROJECT_OF[bare];
  return PROJECT_KEYS.has(bare) ? bare : null;
}

const CONN = process.env.DATABASE_URL || process.env.POSTGRES_URL ||
             process.env.DATABASE_URL_UNPOOLED || process.env.POSTGRES_URL_NON_POOLING;
const sql = CONN ? neon(CONN) : null;
export const hasDb = !!sql;

export const pendingKey = (email) => `pending:${String(email || "").trim().toLowerCase()}`;

let schemaReady = false;
export async function ensureSchema() {
  if (schemaReady || !sql) return;
  await sql`CREATE TABLE IF NOT EXISTS project_access (
    user_id  TEXT PRIMARY KEY,
    email    TEXT,
    name     TEXT,
    projects TEXT[] NOT NULL DEFAULT '{}',
    updated  TIMESTAMPTZ NOT NULL DEFAULT now())`;
  schemaReady = true;
}

// ---- read path (middleware + /api/access) ---------------------------------
let _cache = null, _at = 0;

/** Allowed project keys for a user, or null = every project (the default).
 *
 *  Fails OPEN on a database error: a Neon hiccup must never lock an authenticated lab member
 *  out of the whole site. This is deliberate and safe — the gate that matters (lab membership)
 *  already passed before we get here; this only narrows an already-authenticated user. */
export async function allowedProjectsFor(userId, { force = false } = {}) {
  if (!sql || !userId) return null;
  const now = Date.now();
  if (force || !_cache || now - _at > 60000) {
    try {
      await ensureSchema();
      const rows = await sql`SELECT user_id, projects FROM project_access`;
      const m = {};
      for (const r of rows) m[r.user_id] = r.projects;
      _cache = m; _at = now;
    } catch (_) {
      if (!_cache) _cache = {};        // keep last-known good; empty = everyone sees all
    }
  }
  const v = _cache[userId];
  return Array.isArray(v) ? v : null;
}

export function invalidateCache() { _cache = null; _at = 0; }

// ---- write path -----------------------------------------------------------
/** Called on every successful sign-in. Refreshes the stored display identity, and rebinds an
 *  email pre-seed (`pending:<email>`) onto the real user id the first time we see them.
 *  Best-effort: never throws in a way that could block a legitimate sign-in. */
export async function bindAccessRow(userId, email, name) {
  if (!sql || !userId) return;
  await ensureSchema();
  const existing = await sql`SELECT user_id FROM project_access WHERE user_id = ${userId}`;
  if (existing.length) {
    await sql`UPDATE project_access SET email = ${email || null}, name = ${name || null},
              updated = now() WHERE user_id = ${userId}`;
    invalidateCache();
    return;
  }
  if (!email) return;
  const pk = pendingKey(email);
  const seeded = await sql`SELECT projects FROM project_access WHERE user_id = ${pk}`;
  if (!seeded.length) return;                       // no restriction for this person
  await sql`INSERT INTO project_access (user_id, email, name, projects, updated)
            VALUES (${userId}, ${email}, ${name || null}, ${seeded[0].projects}, now())
            ON CONFLICT (user_id) DO UPDATE SET projects = ${seeded[0].projects}, updated = now()`;
  await sql`DELETE FROM project_access WHERE user_id = ${pk}`;
  invalidateCache();
}

export async function setProjects(key, projects, { email = null, name = null } = {}) {
  if (!sql) throw new Error("no database");
  await ensureSchema();
  const clean = (Array.isArray(projects) ? projects : []).filter((p) => PROJECT_KEYS.has(p));
  await sql`INSERT INTO project_access (user_id, email, name, projects, updated)
            VALUES (${key}, ${email}, ${name}, ${clean}, now())
            ON CONFLICT (user_id) DO UPDATE SET projects = ${clean},
              email = COALESCE(${email}, project_access.email),
              name  = COALESCE(${name},  project_access.name), updated = now()`;
  invalidateCache();
}

/** Remove the restriction entirely — that person goes back to seeing every project. */
export async function clearProjects(key) {
  if (!sql) throw new Error("no database");
  await ensureSchema();
  await sql`DELETE FROM project_access WHERE user_id = ${key}`;
  invalidateCache();
}

export async function allAccessRows() {
  if (!sql) return [];
  await ensureSchema();
  return sql`SELECT user_id, email, name, projects, updated FROM project_access
             ORDER BY COALESCE(name, email, user_id)`;
}
