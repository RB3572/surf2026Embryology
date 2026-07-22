// GET /api/admin — usage analytics for the admin console. Admin-gated twice: the edge
// middleware 404s /api/admin* for non-admins, and this handler re-checks the cookie role
// (defense in depth) returning 404 — never revealing that the endpoint exists.
import { neon } from "@neondatabase/serverless";
import { cookieVal, accountByToken } from "../accounts.mjs";

const CONN = process.env.DATABASE_URL || process.env.POSTGRES_URL ||
             process.env.DATABASE_URL_UNPOOLED || process.env.POSTGRES_URL_NON_POOLING;
const sql = CONN ? neon(CONN) : null;

export default async function handler(req, res) {
  const acct = await accountByToken(cookieVal(req.headers.cookie, "surf_gate"));
  if (!acct || acct.role !== "admin") { res.status(404).json({ error: "Not Found" }); return; }
  if (!sql) {
    res.status(200).json({ ok: false, err: "no-database-url — check the Neon env vars in Vercel",
      users: [], projects: [], recent: [], genes: [], downloads: [] });
    return;
  }
  try {
    // If nothing has been tracked yet the table may not exist; treat that as "empty", not an error.
    const exists = await sql`SELECT to_regclass('public.events') AS t`;
    if (!exists[0] || !exists[0].t) {
      res.status(200).json({ ok: true, empty: true, users: [], projects: [], recent: [], genes: [], downloads: [] });
      return;
    }
    // Optional ?user=<login name> filter. The `users` overview is ALWAYS the full list (it feeds the
    // picker); everything else is scoped to the chosen user. `(${u}::text IS NULL OR usr = ${u})`
    // matches everyone when no filter is set. `u` is parameterized, so it is injection-safe.
    const u = (req.query && typeof req.query.user === "string" && req.query.user.trim()) || null;
    const [users, projects, recent, genes, downloads, totals] = await Promise.all([
      sql`SELECT usr,
            COUNT(*)::int AS events,
            COUNT(*) FILTER (WHERE action='view')::int AS views,
            COUNT(*) FILTER (WHERE action='download')::int AS downloads,
            COUNT(DISTINCT project)::int AS projects,
            COUNT(DISTINCT date_trunc('day', ts))::int AS active_days,
            MAX(ts) AS last_seen, MIN(ts) AS first_seen
          FROM events GROUP BY usr ORDER BY events DESC`,
      sql`SELECT project,
            COUNT(*)::int AS events,
            COUNT(*) FILTER (WHERE action='view')::int AS views,
            COUNT(DISTINCT usr)::int AS users
          FROM events WHERE (${u}::text IS NULL OR usr = ${u})
          GROUP BY project ORDER BY events DESC`,
      sql`SELECT ts, usr, project, action, detail, path FROM events
          WHERE (${u}::text IS NULL OR usr = ${u}) ORDER BY ts DESC LIMIT 800`,
      sql`SELECT detail->>'gene' AS gene, COUNT(*)::int AS n, COUNT(DISTINCT usr)::int AS users
          FROM events WHERE action='gene' AND detail->>'gene' IS NOT NULL AND detail->>'gene' <> ''
            AND (${u}::text IS NULL OR usr = ${u})
          GROUP BY gene ORDER BY n DESC LIMIT 120`,
      sql`SELECT ts, usr, project, detail FROM events
          WHERE action='download' AND (${u}::text IS NULL OR usr = ${u}) ORDER BY ts DESC LIMIT 200`,
      sql`SELECT COUNT(*)::int AS events, COUNT(DISTINCT usr)::int AS users
          FROM events WHERE (${u}::text IS NULL OR usr = ${u})`,
    ]);
    res.status(200).json({ ok: true, filterUser: u, totals: totals[0], users, projects, recent, genes, downloads });
  } catch (e) {
    res.status(200).json({ ok: false, err: String((e && e.message) || e).slice(0, 300),
      users: [], projects: [], recent: [], genes: [], downloads: [] });
  }
}
