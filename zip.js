/* A minimal ZIP writer, STORE method (no compression).
 *
 * The images this packs are already-compressed PNG/JPEG, so DEFLATE would spend time to save
 * nothing. STORE keeps the whole thing to a CRC table and three fixed-layout records, which is
 * far less than pulling a zip library into a page whose only other dependency is Plotly.
 *
 * Loads as a plain script (window.SurfZip) or as an ES module, so the same bytes the browser
 * ships are what the test suite validates.
 */
(function (root, factory) {
  const api = factory();
  if (typeof module === "object" && module.exports) module.exports = api;
  if (root) root.SurfZip = api;
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
  "use strict";

  const TABLE = (() => {
    const t = new Uint32Array(256);
    for (let i = 0; i < 256; i++) {
      let c = i;
      for (let k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1;
      t[i] = c >>> 0;
    }
    return t;
  })();

  // General-purpose bit 11. Names are written as UTF-8; without this flag an unzipper is
  // entitled to read them as CP437, which turns "µm" into mojibake.
  const UTF8_FLAG = 0x0800;

  function crc32(bytes) {
    let c = 0xFFFFFFFF;
    for (let i = 0; i < bytes.length; i++) c = TABLE[(c ^ bytes[i]) & 0xFF] ^ (c >>> 8);
    return (c ^ 0xFFFFFFFF) >>> 0;
  }

  /**
   * @param {{name: string, data: Uint8Array}[]} files
   * @returns {Uint8Array[]} the parts, in order — join them to get the archive
   */
  function zipParts(files) {
    const enc = new TextEncoder(), parts = [], central = [];
    let off = 0;
    for (const f of files) {
      const nm = enc.encode(f.name), crc = crc32(f.data), n = f.data.length;

      const h = new DataView(new ArrayBuffer(30));         // local file header
      h.setUint32(0, 0x04034b50, true);
      h.setUint16(4, 20, true);                            // version needed
      h.setUint16(6, UTF8_FLAG, true);
      h.setUint32(14, crc, true);
      h.setUint32(18, n, true);                            // compressed size == stored size
      h.setUint32(22, n, true);
      h.setUint16(26, nm.length, true);
      parts.push(new Uint8Array(h.buffer), nm, f.data);

      const c = new DataView(new ArrayBuffer(46));         // central directory record
      c.setUint32(0, 0x02014b50, true);
      c.setUint16(4, 20, true);                            // version made by
      c.setUint16(6, 20, true);                            // version needed
      c.setUint16(8, UTF8_FLAG, true);
      c.setUint32(16, crc, true);
      c.setUint32(20, n, true);
      c.setUint32(24, n, true);
      c.setUint16(28, nm.length, true);
      c.setUint32(42, off, true);                          // offset of the local header
      central.push(new Uint8Array(c.buffer), nm);

      off += 30 + nm.length + n;
    }
    const cdSize = central.reduce((a, b) => a + b.length, 0);
    const e = new DataView(new ArrayBuffer(22));           // end of central directory
    e.setUint32(0, 0x06054b50, true);
    e.setUint16(8, files.length, true);
    e.setUint16(10, files.length, true);
    e.setUint32(12, cdSize, true);
    e.setUint32(16, off, true);
    return [...parts, ...central, new Uint8Array(e.buffer)];
  }

  /** One flat Uint8Array — what a test wants. */
  function zipBytes(files) {
    const parts = zipParts(files);
    const out = new Uint8Array(parts.reduce((a, p) => a + p.length, 0));
    let o = 0;
    for (const p of parts) { out.set(p, o); o += p.length; }
    return out;
  }

  /** A Blob — what the browser wants, without concatenating the payloads first. */
  function zipBlob(files) {
    return new Blob(zipParts(files), { type: "application/zip" });
  }

  return { crc32, zipParts, zipBytes, zipBlob };
});
