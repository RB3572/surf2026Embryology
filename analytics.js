/* Lightweight first-party usage tracker. Loaded on every page; sends events to /api/track,
 * which attributes them to the logged-in user via the (HttpOnly) session cookie — the client
 * never sees or sends any identity. Uses global event delegation so pages need no hooks: just
 * the one <script> tag. Fails silently; never blocks or errors the UI. Admin-only analytics
 * are read back through /api/admin. */
(() => {
  const project = (location.pathname.replace(/^\//, "").replace(/\.html$/, "") || "index");

  function send(action, detail) {
    try {
      const body = JSON.stringify({ action, project, path: location.pathname, detail: detail || {} });
      // sendBeacon survives page unloads; fall back to fetch.
      if (navigator.sendBeacon) {
        navigator.sendBeacon("/api/track", new Blob([body], { type: "application/json" }));
      } else {
        fetch("/api/track", { method: "POST", headers: { "content-type": "application/json" },
          body, keepalive: true, credentials: "same-origin" }).catch(() => {});
      }
    } catch (_) { /* never surface */ }
  }
  window.surfTrack = send;                          // pages may call surfTrack('event', {...})

  send("view", { title: document.title });          // page / project view

  // gene navigation — any <select> whose id mentions "gene"
  let lastGene = "";
  document.addEventListener("change", (e) => {
    const t = e.target;
    if (t && t.tagName === "SELECT" && /gene/i.test(t.id || "")) {
      const v = t.value || "";
      if (v && v !== lastGene) { lastGene = v; send("gene", { select: t.id, gene: v }); }
    }
  }, true);

  // chart / CSV downloads — the download & export buttons across the viewers
  document.addEventListener("click", (e) => {
    const b = e.target.closest && e.target.closest(
      ".plot-download, [id*='download'], [id*='csv'], [id*='export'], a[download]");
    if (!b) return;
    const g = document.querySelector("#gene-select, [id^='gene-']");
    send("download", { control: (b.id || b.className || "").slice(0, 60), gene: g ? g.value : null });
  }, true);
})();
