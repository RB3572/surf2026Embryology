// POST /api/track — record one usage event, attributed to the session's user via the
// surf_gate cookie. Open to any authenticated user (that's how Kathy's activity is logged).
// Always returns 200, even on DB failure, so client-side tracking never disrupts the site.
//
// The token→identity map MUST be kept in sync with middleware.js (single source of the
// passwords; here we only need token→{user,role}, which is not secret to reproduce).
import { neon } from "@neondatabase/serverless";

const BY_TOKEN = new Map([
  ["s1_kathytam_b7f3a92c8d1e4056a19d", { user: "Kathy Tam", role: "user" }],
  ["adm_owner_5c1e9a2f7b3d8064c2a7f1", { user: "Admin", role: "admin" }],
]);
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

function cookieVal(header, name) {
  for (const c of (header || "").split(/;\s*/)) {
    const i = c.indexOf("=");
    if (i > 0 && c.slice(0, i) === name) return c.slice(i + 1);
  }
  return null;
}

export default async function handler(req, res) {
  if (req.method !== "POST") { res.status(405).json({ ok: false }); return; }
  const token = cookieVal(req.headers.cookie, "surf_gate");
  const acct = token && BY_TOKEN.get(token);
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
