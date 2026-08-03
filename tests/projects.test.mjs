// Tests for the per-project access helpers. These are pure (no DB) — the database-backed
// functions in lib/projects.mjs are exercised against Neon in the deployed environment.

import test from "node:test";
import assert from "node:assert/strict";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

import { PROJECTS, PROJECT_KEYS, PROJECT_OF, projectForPath, pendingKey } from "../lib/projects.mjs";

const ROOT = path.join(path.dirname(fileURLToPath(import.meta.url)), "..");

test("every project key maps to a page that actually exists", () => {
  for (const p of PROJECTS) {
    assert.ok(fs.existsSync(path.join(ROOT, `${p.key}.html`)),
      `${p.key}.html is listed in PROJECTS but missing from the repo`);
  }
});

test("every project page on disk is listed in PROJECTS", () => {
  // A page that isn't listed can never be restricted — it would be silently visible to
  // everyone, which is exactly the drift the old hand-maintained list suffered from.
  const skip = new Set(["index.html", "admin.html", "admin.card.html", "login.html"]);
  const pages = fs.readdirSync(ROOT)
    .filter((f) => f.endsWith(".html") && !skip.has(f));
  for (const f of pages) {
    assert.ok(PROJECT_KEYS.has(f.replace(/\.html$/, "")),
      `${f} exists but is not listed in PROJECTS (it could not be access-controlled)`);
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
