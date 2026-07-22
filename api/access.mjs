// GET/POST /api/access — per-user project access.
//   GET               (any authed user)  → { me: [projects] | null }   (null = all projects)
//   GET ?matrix=1     (admin only)        → { projects:[{key,label}], users:[{usr,role,projects|null}] }
//   POST {usr,projects}(admin only)       → upsert that user's allowed-project list ([] = none, absent row = all)
//
// Stored in the Neon `access` table (usr TEXT PK, projects TEXT[]). A user with NO row sees every
// project; a row restricts them to exactly that list. Kathy is seeded to the two pronuclei projects.
// Identity (base + admin-added logins) comes from the shared accounts module.
import { neon } from "@neondatabase/serverless";
import { cookieVal, accountByToken, allAccounts } from "../accounts.mjs";

// The canonical project list (key = the page's basename without .html). Keep in sync with the
// landing cards + middleware PROJECT_OF.
export const PROJECTS = [
  { key: "pronuclei", label: "Pronuclei Distance vs Transcripts" },
  { key: "extpt", label: "Extended Pseudotime Analysis" },
  { key: "segments", label: "Segment Gene Enrichment" },
  { key: "zygote-planes", label: "Zygote Division Planes" },
  { key: "sperm-division", label: "Sperm Division Plane" },
  { key: "pronuclei-assignments", label: "Pronuclei Assignments" },
  { key: "diffusion", label: "Gene Diffusion Rates" },
  { key: "sperm-pca", label: "Sperm Prediction (PCA)" },
  { key: "axes", label: "Fertilization Geometry" },
  { key: "alphabeta", label: "Sperm α/β · 2-cell" },
];
const KEYS = new Set(PROJECTS.map((p) => p.key));
const DEFAULTS = { "Kathy Tam": ["pronuclei", "extpt"] };   // seeded on first run

const CONN = process.env.DATABASE_URL || process.env.POSTGRES_URL ||
             process.env.DATABASE_URL_UNPOOLED || process.env.POSTGRES_URL_NON_POOLING;
const sql = CONN ? neon(CONN) : null;
let ready = false;

async function ensureSchema() {
  if (ready || !sql) return;
  await sql`CREATE TABLE IF NOT EXISTS access (usr TEXT PRIMARY KEY, projects TEXT[])`;
  const n = await sql`SELECT COUNT(*)::int AS c FROM access`;
  if (!n[0].c) for (const [u, ps] of Object.entries(DEFAULTS))
    await sql`INSERT INTO access (usr, projects) VALUES (${u}, ${ps}) ON CONFLICT (usr) DO NOTHING`;
  ready = true;
}
export default async function handler(req, res) {
  const acct = await accountByToken(cookieVal(req.headers.cookie, "surf_gate"));
  if (!acct) { res.status(401).json({ ok: false }); return; }
  if (!sql) { res.status(200).json({ ok: false, err: "no-database-url", me: null, projects: PROJECTS }); return; }
  try {
    await ensureSchema();
    if (req.method === "POST") {
      if (acct.role !== "admin") { res.status(404).json({ error: "Not Found" }); return; }
      let b = req.body; if (typeof b === "string") { try { b = JSON.parse(b); } catch { b = {}; } }
      const usr = String((b && b.usr) || "");
      const projects = Array.isArray(b && b.projects) ? b.projects.filter((p) => KEYS.has(p)) : [];
      const accs = await allAccounts();
      if (!usr || !accs.some((u) => u.user === usr)) { res.status(400).json({ ok: false, err: "unknown user" }); return; }
      await sql`INSERT INTO access (usr, projects) VALUES (${usr}, ${projects})
                ON CONFLICT (usr) DO UPDATE SET projects = ${projects}`;
      res.status(200).json({ ok: true }); return;
    }
    const rows = await sql`SELECT usr, projects FROM access`;
    const map = {}; rows.forEach((r) => (map[r.usr] = r.projects));
    if (req.query && req.query.matrix) {
      if (acct.role !== "admin") { res.status(404).json({ error: "Not Found" }); return; }
      const accs = await allAccounts();
      res.status(200).json({ ok: true, projects: PROJECTS,
        users: accs.map((u) => ({ usr: u.user, role: u.role, projects: map[u.user] == null ? null : map[u.user] })) });
      return;
    }
    // default: the caller's own access (null = all projects)
    const mine = map[acct.user];
    res.status(200).json({ ok: true, me: mine == null ? null : mine, admin: acct.role === "admin" });
  } catch (e) {
    res.status(200).json({ ok: false, err: String((e && e.message) || e).slice(0, 200), me: null, projects: PROJECTS });
  }
}
