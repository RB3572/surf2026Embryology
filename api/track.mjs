// POST /api/track — record one usage event, attributed to the session's user via the
// surf_gate cookie. Open to any authenticated user (that's how Kathy's activity is logged).
// Always returns 200, even on DB failure, so client-side tracking never disrupts the site.
//
// Identity comes from the shared accounts module (same source as the gate) — no local copy to drift.
import { neon } from "@neondatabase/serverless";
import { cookieVal, accountByToken } from "../accounts.mjs";

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
  await sql`CREATE INDEX IF NOT EXISTS events_ts_idx ON events (ts DESC)`;
  await sql`CREATE INDEX IF NOT EXISTS events_usr_idx ON events (usr)`;
  schemaReady = true;
}

export default async function handler(req, res) {
  if (req.method !== "POST") { res.status(405).json({ ok: false }); return; }
  const token = cookieVal(req.headers.cookie, "surf_gate");
  const acct = await accountByToken(token);
  if (!acct) { res.status(401).json({ ok: false }); return; }
  if (!sql) { res.status(200).json({ ok: false, err: "no-database-url" }); return; }
  try {
    await ensureSchema();
    let b = req.body;
    if (typeof b === "string") { try { b = JSON.parse(b); } catch { b = {}; } }
    b = b || {};
    await sql`INSERT INTO events (token, usr, role, project, action, detail, path, ua)
      VALUES (${token}, ${acct.user}, ${acct.role},
              ${String(b.project || "").slice(0, 80)}, ${String(b.action || "").slice(0, 40)},
              ${JSON.stringify(b.detail || {})}, ${String(b.path || "").slice(0, 200)},
              ${String(req.headers["user-agent"] || "").slice(0, 200)})`;
    res.status(200).json({ ok: true });
  } catch (e) {
    res.status(200).json({ ok: false, err: String((e && e.message) || e).slice(0, 200) });
  }
}
