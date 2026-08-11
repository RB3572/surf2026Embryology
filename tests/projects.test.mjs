// Tests for the per-project access helpers. These are pure (no DB) — the database-backed
// functions in lib/projects.mjs are exercised against Neon in the deployed environment.

import test from "node:test";
import assert from "node:assert/strict";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

import { PROJECTS, PROJECT_KEYS, PROJECT_OF, projectForPath, pendingKey, ADMIN_ONLY_PAGES } from "../lib/projects.mjs";

const ROOT = path.join(path.dirname(fileURLToPath(import.meta.url)), "..");

test("every project key maps to a page that actually exists", () => {
  for (const p of PROJECTS) {
    assert.ok(fs.existsSync(path.join(ROOT, `${p.key}.html`)),
      `${p.key}.html is listed in PROJECTS but missing from the repo`);
  }
});

test("every project page on disk is access-controlled", () => {
  // A page in NEITHER list can never be restricted — it would be silently visible to everyone,
  // which is exactly the drift the old hand-maintained list suffered from. There are two tiers:
  // PROJECTS is the per-user allowlist, ADMIN_ONLY_PAGES is the stricter admin-only tier that is
  // deliberately kept out of PROJECTS so the access matrix doesn't offer a pointless toggle.
  // Either tier is access control; being in neither is the hole this guards.
  const skip = new Set(["index.html", "admin.html", "admin.card.html", "login.html"]);
  const pages = fs.readdirSync(ROOT)
    .filter((f) => f.endsWith(".html") && !skip.has(f));
  for (const f of pages) {
    const key = f.replace(/\.html$/, "");
    assert.ok(PROJECT_KEYS.has(key) || ADMIN_ONLY_PAGES.has(key),
      `${f} exists but is in neither PROJECTS nor ADMIN_ONLY_PAGES (it could not be access-controlled)`);
  }
});

test("the two access tiers are disjoint", () => {
  // A key in both would be ambiguous: the middleware would have to choose whether the per-user
  // allowlist or the admin-only rule wins, and the admin console would show a toggle that does
  // nothing. Keeping them disjoint means the tier a page belongs to is never in question.
  const both = [...ADMIN_ONLY_PAGES].filter((k) => PROJECT_KEYS.has(k));
  assert.deepEqual(both, [], `keys in both tiers: ${both.join(", ")}`);
});

test("every admin-only page exists on disk", () => {
  for (const key of ADMIN_ONLY_PAGES) {
    assert.ok(fs.existsSync(path.join(ROOT, `${key}.html`)),
      `${key}.html is listed in ADMIN_ONLY_PAGES but missing from the repo`);
  }
});

test("project keys are unique", () => {
  assert.equal(PROJECT_KEYS.size, PROJECTS.length);
});

test("PROJECT_OF is derived from PROJECTS and cannot drift", () => {
  assert.equal(Object.keys(PROJECT_OF).length, PROJECTS.length);
  for (const p of PROJECTS) assert.equal(PROJECT_OF[`${p.key}.html`], p.key);
});

test("projectForPath matches both the .html and extensionless forms", () => {
  assert.equal(projectForPath("/sperm-sphere.html"), "sperm-sphere");
  assert.equal(projectForPath("/sperm-sphere"), "sperm-sphere");
  assert.equal(projectForPath("sperm-sphere.html"), "sperm-sphere");
});

test("projectForPath returns null for non-project paths", () => {
  for (const p of ["/", "/index.html", "/admin", "/data/zygote_sperm.json",
                   "/viewer-core.js", "/nope.html", "", null, undefined]) {
    assert.equal(projectForPath(p), null, `${p} should not be treated as a project`);
  }
});

test("pendingKey normalises the email so a pre-seed always rebinds", () => {
  assert.equal(pendingKey("A.Person@Example.COM "), "pending:a.person@example.com");
  assert.equal(pendingKey("  x@y.z"), "pending:x@y.z");
});
