// GET /api/admin — usage analytics for the admin console. Admin-gated twice: the edge
// middleware 404s /api/admin* for non-admins, and this handler re-checks the session's admin
// flag (defense in depth) returning 404 — never revealing that the endpoint exists.
//
// Admin = listed in ADMIN_USER_IDS, or a lab admin/PI in Lab Logger (see lib/session.mjs
// isAdminUser). There are no passwords and no local roles any more.
//
// People are keyed on COALESCE(user_id, usr): rows written since the SSO migration carry the
// Supabase user id, while pre-SSO rows only have a display name. Both are kept so the history
// survives; a person who used the site before AND after may appear under both keys.
import { neon } from "@neondatabase/serverless";
import { SESSION_COOKIE, verifyToken, readCookie } from "../lib/session.mjs";

const CONN = process.env.DATABASE_URL || process.env.POSTGRES_URL ||
             process.env.DATABASE_URL_UNPOOLED || process.env.POSTGRES_URL_NON_POOLING;
const sql = CONN ? neon(CONN) : null;

export default async function handler(req, res) {
  const session = await verifyToken(
    readCookie(req.headers.cookie, SESSION_COOKIE), process.env.SESSION_SECRET);

  res.setHeader("content-type", "application/json");
  res.setHeader("Cache-Control", "no-store");

  if (!session || !session.adm) {
    res.statusCode = 404;
    res.end('{"error":"Not Found"}');
    return;
  }

  const empty = { ok: true, empty: true, viewer: { name: session.nm, email: session.em },
                  users: [], projects: [], recent: [], genes: [], downloads: [] };

  if (!sql) {
    res.statusCode = 200;
    res.end(JSON.stringify({ ...empty, ok: false,
      err: "no-database-url — check the Neon env vars in Vercel" }));
    return;
  }

  try {
    // If nothing has been tracked yet the table may not exist; treat that as "empty".
    const exists = await sql`SELECT to_regclass('public.events') AS t`;
    if (!exists[0] || !exists[0].t) {
      res.statusCode = 200; res.end(JSON.stringify(empty)); return;
    }
    // Make sure the SSO columns exist even if no event has been written since the migration
    // (otherwise the queries below would fail on a pre-migration table).
    await sql`ALTER TABLE events ADD COLUMN IF NOT EXISTS user_id TEXT`;
    await sql`ALTER TABLE events ADD COLUMN IF NOT EXISTS email TEXT`;

    // Optional ?user=<person key> filter. `u` is parameterized, so it is injection-safe.
    const u = (req.query && typeof req.query.user === "string" && req.query.user.trim()) || null;
    const [users, projects, recent, genes, downloads, totals] = await Promise.all([
      // Always the full roster — this feeds the picker.
      sql`SELECT COALESCE(user_id, usr) AS person, MAX(usr) AS usr, MAX(email) AS email,
            COUNT(*)::int AS events,
            COUNT(*) FILTER (WHERE action='view')::int AS views,
            COUNT(*) FILTER (WHERE action='download')::int AS downloads,
            COUNT(DISTINCT project)::int AS projects,
            COUNT(DISTINCT date_trunc('day', ts))::int AS active_days,
            MAX(ts) AS last_seen, MIN(ts) AS first_seen
          FROM events GROUP BY COALESCE(user_id, usr) ORDER BY events DESC`,
      sql`SELECT project,
            COUNT(*)::int AS events,
            COUNT(*) FILTER (WHERE action='view')::int AS views,
            COUNT(DISTINCT COALESCE(user_id, usr))::int AS users
          FROM events WHERE (${u}::text IS NULL OR COALESCE(user_id, usr) = ${u})
          GROUP BY project ORDER BY events DESC`,
      sql`SELECT ts, usr, project, action, detail, path FROM events
          WHERE (${u}::text IS NULL OR COALESCE(user_id, usr) = ${u})
          ORDER BY ts DESC LIMIT 800`,
      sql`SELECT detail->>'gene' AS gene, COUNT(*)::int AS n,
            COUNT(DISTINCT COALESCE(user_id, usr))::int AS users
          FROM events WHERE action='gene' AND detail->>'gene' IS NOT NULL AND detail->>'gene' <> ''
            AND (${u}::text IS NULL OR COALESCE(user_id, usr) = ${u})
          GROUP BY gene ORDER BY n DESC LIMIT 120`,
      sql`SELECT ts, usr, project, detail FROM events
          WHERE action='download' AND (${u}::text IS NULL OR COALESCE(user_id, usr) = ${u})
          ORDER BY ts DESC LIMIT 200`,
      sql`SELECT COUNT(*)::int AS events, COUNT(DISTINCT COALESCE(user_id, usr))::int AS users
          FROM events WHERE (${u}::text IS NULL OR COALESCE(user_id, usr) = ${u})`,
    ]);
    res.statusCode = 200;
    res.end(JSON.stringify({ ok: true, filterUser: u, viewer: { name: session.nm, email: session.em },
      totals: totals[0], users, projects, recent, genes, downloads }));
  } catch (e) {
    res.statusCode = 200;
    res.end(JSON.stringify({ ...empty, ok: false,
      err: String((e && e.message) || e).slice(0, 300) }));
  }
}
