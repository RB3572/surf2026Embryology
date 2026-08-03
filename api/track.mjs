// POST /api/track — record one usage event, attributed to the signed-in user via the httpOnly
// session cookie. The client never sends identity; it can only say WHAT happened, never WHO.
// Always returns 200 (even on DB failure) so tracking can never break the site.
//
// Identity is now the Supabase user id from Lab Logger SSO. The pre-SSO rows in this table are
// keyed only by display name (`usr`), so both columns are written and the admin console keys
// people on COALESCE(user_id, usr) — old history is preserved rather than orphaned.
import { neon } from "@neondatabase/serverless";
import { SESSION_COOKIE, verifyToken, readCookie } from "../lib/session.mjs";

const CONN = process.env.DATABASE_URL || process.env.POSTGRES_URL ||
             process.env.DATABASE_URL_UNPOOLED || process.env.POSTGRES_URL_NON_POOLING;
const sql = CONN ? neon(CONN) : null;
let schemaReady = false;

async function ensureSchema() {
  if (schemaReady || !sql) return;
  await sql`CREATE TABLE IF NOT EXISTS events (
    id BIGSERIAL PRIMARY KEY,
    ts TIMESTAMPTZ NOT NULL DEFAULT now(),
    token TEXT, usr TEXT, role TEXT,
    project TEXT, action TEXT, detail JSONB, path TEXT, ua TEXT)`;
  // Added by the SSO migration; the pre-existing table won't have these.
  await sql`ALTER TABLE events ADD COLUMN IF NOT EXISTS user_id TEXT`;
  await sql`ALTER TABLE events ADD COLUMN IF NOT EXISTS email TEXT`;
  await sql`CREATE INDEX IF NOT EXISTS events_ts_idx ON events (ts DESC)`;
  await sql`CREATE INDEX IF NOT EXISTS events_usr_idx ON events (usr)`;
  await sql`CREATE INDEX IF NOT EXISTS events_user_id_idx ON events (user_id)`;
  schemaReady = true;
}

const clip = (v, n) => (v == null ? null : String(v).slice(0, n));

export default async function handler(req, res) {
  res.setHeader("content-type", "application/json");
  res.setHeader("Cache-Control", "no-store");

  if (req.method !== "POST") { res.statusCode = 405; res.end('{"ok":false}'); return; }

  const session = await verifyToken(
    readCookie(req.headers.cookie, SESSION_COOKIE), process.env.SESSION_SECRET);
  if (!session) { res.statusCode = 401; res.end('{"ok":false}'); return; }

  if (!sql) { res.statusCode = 200; res.end('{"ok":false,"err":"no-database-url"}'); return; }

  try {
    await ensureSchema();
    let b = req.body;
    if (typeof b === "string") { try { b = JSON.parse(b); } catch { b = {}; } }
    b = b || {};
    await sql`INSERT INTO events (user_id, email, usr, role, project, action, detail, path, ua)
      VALUES (${session.sub}, ${clip(session.em, 200)}, ${clip(session.nm, 120)},
              ${session.adm ? "admin" : (clip(session.rl, 20) || "user")},
              ${clip(b.project, 80)}, ${clip(b.action, 40)},
              ${JSON.stringify(b.detail || {})}, ${clip(b.path, 200)},
              ${clip(req.headers["user-agent"], 200)})`;
    res.statusCode = 200;
    res.end('{"ok":true}');
  } catch (e) {
    res.statusCode = 200;      // never disrupt the site
    res.end(JSON.stringify({ ok: false, err: String((e && e.message) || e).slice(0, 200) }));
  }
}
