/* zip.js writes a binary container by hand, so "it downloaded something" is not evidence it
 * works. These build archives with the shipped code and hand the bytes to Python's zipfile,
 * which is an independent implementation — if the headers, offsets or CRCs are wrong it fails
 * to open or the payloads come back changed.
 */
import { execFileSync } from "node:child_process";
import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";
import assert from "node:assert/strict";
import path from "node:path";
import test from "node:test";

const HERE = path.dirname(fileURLToPath(import.meta.url));
// The package is "type": "module", so Node parses zip.js as ESM; the UMD wrapper then takes its
// global branch rather than the CommonJS one. Load it, then read what it published.
const require = createRequire(import.meta.url);
const loaded = require(path.join(HERE, "..", "zip.js"));
const { crc32, zipBytes } = globalThis.SurfZip || loaded;

const enc = (s) => new TextEncoder().encode(s);

/** Round-trip through Python: returns {names, contents, bad} as reported by zipfile. */
function readBackWithPython(bytes) {
  const py = `
import base64, io, json, sys, zipfile
raw = base64.b64decode(sys.stdin.read())
z = zipfile.ZipFile(io.BytesIO(raw))
bad = z.testzip()
print(json.dumps({
    "names": z.namelist(),
    "contents": {n: base64.b64encode(z.read(n)).decode() for n in z.namelist()},
    "bad": bad,
}))
`;
  const out = execFileSync("python3", ["-c", py], {
    input: Buffer.from(bytes).toString("base64"),
    encoding: "utf8",
    maxBuffer: 64 * 1024 * 1024,
  });
  return JSON.parse(out);
}

test("crc32 matches the known CRC-32 of '123456789'", () => {
  assert.equal(crc32(enc("123456789")), 0xCBF43926);
});

test("crc32 of the empty input is 0", () => {
  assert.equal(crc32(new Uint8Array(0)), 0);
});

test("python zipfile opens the archive and the payloads survive byte for byte", () => {
  const files = [
    { name: "README.txt", data: enc("embryo\tstage\nZ-P3-fov23\tzygote\n") },
    { name: "a_488_zslice.png", data: new Uint8Array([137, 80, 78, 71, 13, 10, 26, 10, 0, 1, 2, 255]) },
    { name: "nested/dir/405_maxZ.jpg", data: new Uint8Array(3000).map((_, i) => (i * 7) & 0xFF) },
  ];
  const got = readBackWithPython(zipBytes(files));

  assert.equal(got.bad, null, "zipfile.testzip() found a corrupt member");
  assert.deepEqual(got.names, files.map((f) => f.name));
  for (const f of files) {
    const back = Buffer.from(got.contents[f.name], "base64");
    assert.deepEqual(new Uint8Array(back), f.data, `payload changed for ${f.name}`);
  }
});

test("an empty archive is still a valid, openable zip", () => {
  const got = readBackWithPython(zipBytes([]));
  assert.deepEqual(got.names, []);
});

test("non-ASCII names survive (embryo ids and µm show up in filenames)", () => {
  const files = [{ name: "distance_µm_Z–P3.txt", data: enc("2.4 µm") }];
  const got = readBackWithPython(zipBytes(files));
  assert.deepEqual(got.names, ["distance_µm_Z–P3.txt"]);
  assert.equal(Buffer.from(got.contents[got.names[0]], "base64").toString("utf8"), "2.4 µm");
});

test("local header offsets stay correct as entries accumulate", () => {
  // 40 entries of differing length: a mistake in the running offset shows up as a failure to
  // open the later members rather than the earlier ones.
  const files = Array.from({ length: 40 }, (_, i) => ({
    name: `f${i}.bin`,
    data: new Uint8Array(i * 137).fill(i & 0xFF),
  }));
  const got = readBackWithPython(zipBytes(files));
  assert.equal(got.bad, null);
  assert.equal(got.names.length, 40);
  assert.equal(Buffer.from(got.contents["f39.bin"], "base64").length, 39 * 137);
});
