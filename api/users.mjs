// Admin-only login management for the console.
//   GET            → { accounts: [{usr, role, pw, kind}] }  (all logins WITH passwords, so the
//                    admin can read them back to hand out)
//   POST {usr,pw,role}  → add (or update) a login in the Neon `users` table
//   DELETE ?usr=… (or {usr})  → remove an added login (built-in ones can't be removed)
// Double-gated: the edge middleware 404s /api/* admin surfaces for non-admins, and this re-checks
// the cookie role and returns 404 (never reveals the endpoint exists).
import { cookieVal, accountByToken, BASE, dbUsers, addUser, removeUser, hasDb } from "../accounts.mjs";

export default async function handler(req, res) {
  const acct = await accountByToken(cookieVal(req.headers.cookie, "surf_gate"));
  if (!acct || acct.role !== "admin") { res.status(404).json({ error: "Not Found" }); return; }
  try {
    if (req.method === "POST") {
      let b = req.body; if (typeof b === "string") { try { b = JSON.parse(b); } catch { b = {}; } }
      b = b || {};
      const usr = String(b.usr || "").trim();
      const pw = String(b.pw || "");
      const role = b.role === "admin" ? "admin" : "user";
      if (!usr || !pw) { res.status(400).json({ ok: false, err: "login name and password are both required" }); return; }
      if (BASE.some((a) => a.user.toLowerCase() === usr.toLowerCase())) {
        res.status(400).json({ ok: false, err: `"${usr}" is a built-in login — pick another name` }); return; }
      if (!hasDb) { res.status(200).json({ ok: false, err: "no database configured — can't add logins" }); return; }
      await addUser(usr, pw, role);
      res.status(200).json({ ok: true }); return;
    }
    if (req.method === "DELETE") {
      let usr = (req.query && req.query.usr) || "";
      if (!usr) { let b = req.body; if (typeof b === "string") { try { b = JSON.parse(b); } catch { b = {}; } } usr = (b && b.usr) || ""; }
      usr = String(usr).trim();
      if (!usr) { res.status(400).json({ ok: false, err: "which login?" }); return; }
      if (BASE.some((a) => a.user === usr)) { res.status(400).json({ ok: false, err: "built-in logins can't be removed" }); return; }
      if (!hasDb) { res.status(200).json({ ok: false, err: "no database configured" }); return; }
      await removeUser(usr);
      res.status(200).json({ ok: true }); return;
    }
    // GET — every login with its password (admin-only)
    const db = await dbUsers(true);
    const accounts = BASE.map((a) => ({ usr: a.user, role: a.role, pw: a.pw, kind: "built-in" }))
      .concat(db.map((a) => ({ usr: a.user, role: a.role, pw: a.pw, kind: "added" })));
    res.status(200).json({ ok: true, accounts, hasDb });
  } catch (e) {
    res.status(200).json({ ok: false, err: String((e && e.message) || e).slice(0, 200) });
  }
}
