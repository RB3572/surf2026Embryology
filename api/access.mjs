// GET/POST/DELETE /api/access — per-person project access.
//
//   GET                (any lab member) → { me: [projects] | null, admin }   null = all projects
//   GET ?matrix=1      (admin only)     → { projects, users:[{key,user_id,email,name,projects,pending}] }
//   POST {key|email, projects}(admin)   → restrict that person to exactly `projects`
//   DELETE ?key=…      (admin only)     → drop the restriction (back to every project)
//
// Keyed on the Supabase user id. An admin may also restrict someone by EMAIL before they have
// ever signed in; that row is rebound to their real user id on first sign-in (see
// lib/projects.mjs bindAccessRow), so the restriction applies from their very first visit.

import { neon } from "@neondatabase/serverless";
import { SESSION_COOKIE, verifyToken, readCookie } from "../lib/session.mjs";
import {
  PROJECTS, PROJECT_KEYS, hasDb, pendingKey,
  allowedProjectsFor, setProjects, clearProjects, allAccessRows,
} from "../lib/projects.mjs";

const CONN = process.env.DATABASE_URL || process.env.POSTGRES_URL ||
             process.env.DATABASE_URL_UNPOOLED || process.env.POSTGRES_URL_NON_POOLING;
const sql = CONN ? neon(CONN) : null;

const json = (res, status, body) => {
  res.statusCode = status;
  res.setHeader("content-type", "application/json");
  res.setHeader("Cache-Control", "no-store");
  res.end(JSON.stringify(body));
};

export default async function handler(req, res) {
  const session = await verifyToken(
    readCookie(req.headers.cookie, SESSION_COOKIE), process.env.SESSION_SECRET);
  if (!session) return json(res, 401, { ok: false });

  if (!hasDb) {
    return json(res, 200, { ok: false, err: "no-database-url", me: null,
                            admin: !!session.adm, projects: PROJECTS });
  }

  try {
    if (req.method === "POST") {
      if (!session.adm) return json(res, 404, { error: "Not Found" });
      let b = req.body;
      if (typeof b === "string") { try { b = JSON.parse(b); } catch { b = {}; } }
      b = b || {};

      // Either an existing key (user id / pending:… ) or a raw email to pre-seed.
      const email = String(b.email || "").trim();
      const key = String(b.key || b.user_id || "").trim() || (email ? pendingKey(email) : "");
      if (!key) return json(res, 400, { ok: false, err: "need a user id or an email" });

      const projects = Array.isArray(b.projects)
        ? b.projects.filter((p) => PROJECT_KEYS.has(p)) : [];
      await setProjects(key, projects, { email: email || null, name: b.name || null });
      return json(res, 200, { ok: true });
    }

    if (req.method === "DELETE") {
      if (!session.adm) return json(res, 404, { error: "Not Found" });
      let key = (req.query && req.query.key) || "";
      if (!key) {
        let b = req.body;
        if (typeof b === "string") { try { b = JSON.parse(b); } catch { b = {}; } }
        key = (b && (b.key || b.user_id)) || "";
      }
      key = String(key).trim();
      if (!key) return json(res, 400, { ok: false, err: "which person?" });
      await clearProjects(key);
      return json(res, 200, { ok: true });
    }

    // ---- GET ?matrix=1 — the admin console's editing grid ----
    if (req.query && req.query.matrix) {
      if (!session.adm) return json(res, 404, { error: "Not Found" });
      const rows = await allAccessRows();
      const byKey = new Map(rows.map((r) => [r.user_id, r]));

      // Everyone we know about = people with a restriction + everyone who has ever visited
      // (the analytics table is where signed-in identities show up).
      let seen = [];
      try {
        const t = await sql`SELECT to_regclass('public.events') AS t`;
        if (t[0]?.t) {
          seen = await sql`SELECT user_id, MAX(email) AS email, MAX(usr) AS name,
                                  MAX(ts) AS last_seen
                           FROM events WHERE user_id IS NOT NULL
                           GROUP BY user_id ORDER BY MAX(ts) DESC`;
        }
      } catch (_) { /* analytics is optional */ }

      const users = [];
      for (const s of seen) {
        const r = byKey.get(s.user_id);
        users.push({
          key: s.user_id, user_id: s.user_id,
          email: r?.email || s.email || "", name: r?.name || s.name || "",
          projects: r ? r.projects : null, pending: false, last_seen: s.last_seen,
        });
        byKey.delete(s.user_id);
      }
      // Rows for people who haven't signed in yet (email pre-seeds), plus any restriction on
      // someone with no analytics rows.
      for (const r of byKey.values()) {
        const isPending = r.user_id.startsWith("pending:");
        users.push({
          key: r.user_id, user_id: isPending ? "" : r.user_id,
          email: r.email || (isPending ? r.user_id.slice("pending:".length) : ""),
          name: r.name || "", projects: r.projects,
          pending: isPending, last_seen: null,
        });
      }
      return json(res, 200, { ok: true, projects: PROJECTS, users });
    }

    // ---- GET — the caller's own access (null = every project) ----
    const mine = session.adm ? null : await allowedProjectsFor(session.sub, { force: true });
    return json(res, 200, { ok: true, me: mine, admin: !!session.adm });
  } catch (e) {
    return json(res, 200, { ok: false, err: String((e && e.message) || e).slice(0, 200),
                            me: null, admin: !!session.adm, projects: PROJECTS });
  }
}
