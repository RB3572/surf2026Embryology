// Vercel Edge Middleware — the real access gate for the SURF 2026 embryology site.
//
// Authorization is delegated entirely to Lab Logger: you are in if you are a member of an
// allowed lab (checked at sign-in against Supabase under the user's own RLS). There is no
// user list and no password in this repo.
//
// This runs BEFORE any file is served, which is the whole point: the valuable part of this
// site is the build_*.py output under /data plus the per-project pages, and those are plain
// static assets. Client-side JS could not protect them; this does.
//
// Unauthenticated:
//   • HTML navigation  -> redirect to /login?next=…
//   • anything else    -> 401 (so fetch()/curl for /data gets a clean refusal)

import { next } from "@vercel/edge";
import {
  SESSION_COOKIE, PROVIDER_COOKIE, AUTOTRY_COOKIE,
  verifyToken, readCookie, normalizeProvider, safeNext,
} from "./lib/session.mjs";
import { projectForPath, allowedProjectsFor, ADMIN_ONLY_PAGES } from "./lib/projects.mjs";

export const config = {
  // Guard EVERY path that serves content — pages, /data, scripts, styles, images — not just
  // HTML. Excluded: the auth endpoints themselves (or sign-in could never complete) and
  // Vercel's internals. /login IS matched, so a direct or bookmarked visit can still
  // auto-continue; it is never redirected to itself (see the guard below).
  matcher: ["/((?!api/auth|_vercel|favicon\\.ico|robots\\.txt).*)"],
};

export default async function middleware(request) {
  const secret = process.env.SESSION_SECRET;

  // FAIL CLOSED: a misconfigured deployment denies rather than serving the dataset.
  if (!secret) {
    return new Response("Auth is not configured on this deployment.", {
      status: 503,
      headers: { "Cache-Control": "no-store", "Content-Type": "text/plain; charset=utf-8" },
    });
  }

  const url = new URL(request.url);
  const session = await verifyToken(
    readCookie(request.headers.get("cookie"), SESSION_COOKIE), secret);

  if (session) {
    const p = url.pathname;

    // The admin console, its API, and a handful of admin-only project pages (see
    // ADMIN_ONLY_PAGES) are all admin-only. Answer 404 rather than 403 so a non-admin member
    // cannot even tell the page exists.
    const bareKey = p.replace(/^\/+/, "").replace(/\/+$/, "").replace(/\.html?$/, "");
    const adminOnly = p === "/admin" || p.startsWith("/admin/") || p.startsWith("/admin.") ||
                      p.startsWith("/api/admin") || ADMIN_ONLY_PAGES.has(bareKey);
    if (adminOnly && !session.adm) {
      return new Response("Not Found", {
        status: 404,
        headers: { "Cache-Control": "no-store", "Content-Type": "text/plain; charset=utf-8" },
      });
    }

    // Per-person project access. Admins always see everything; a restricted member asking for
    // a project outside their list is bounced to the landing page. Fails OPEN on a DB error.
    if (!session.adm) {
      const proj = projectForPath(p);
      if (proj) {
        try {
          const allowed = await allowedProjectsFor(session.sub);
          if (allowed && !allowed.includes(proj)) {
            const to = new URL("/", url);
            to.searchParams.set("denied", proj);
            return Response.redirect(to.toString(), 303);
          }
        } catch (_) { /* fail open: serve */ }
      }
    }

    // Authenticated: serve, but make sure no shared cache ever stores an authenticated
    // response and hands it to somebody else.
    const res = next();
    res.headers.set("Vary", "Cookie");
    return res;
  }

  const wantsHtml = (request.headers.get("accept") || "").includes("text/html");

  if (!wantsHtml) {
    return new Response("Unauthorized", {
      status: 401,
      headers: { "Cache-Control": "no-store", "Content-Type": "text/plain; charset=utf-8" },
    });
  }

  // ---- Continuous sign-in --------------------------------------------------
  // We can't read Lab Logger's Supabase session (different origin, and third-party cookies
  // are dead), so we must bounce through the identity provider. That bounce IS silent when
  // the user already has a Google/Apple session — the only friction is our own provider
  // chooser. So skip the chooser whenever the provider is known: from the incoming link
  // (?provider=google, which Lab Logger can add) or from the last successful sign-in on this
  // browser. Result: arriving from Lab Logger costs zero clicks.
  //
  // We never GUESS. With no signal we show the chooser, because sending an Apple user to
  // Google can silently authenticate a DIFFERENT account that isn't in the lab.
  const cookies = request.headers.get("cookie");
  const hinted = normalizeProvider(url.searchParams.get("provider"));
  const remembered = normalizeProvider(readCookie(cookies, PROVIDER_COOKIE));
  const provider = hinted || remembered;
  // If a silent attempt just failed we'd loop; fall back to the chooser once.
  const alreadyTried = readCookie(cookies, AUTOTRY_COOKIE) === "1";
  // After an explicit sign-out, never bounce them straight back in — otherwise signing out
  // would be impossible.
  const justSignedOut = url.searchParams.has("signedout");
  const onLoginPage = url.pathname === "/login" || url.pathname === "/login.html";

  const nextPath = onLoginPage
    ? safeNext(url.searchParams.get("next"))
    : url.pathname + url.search;

  if (provider && !alreadyTried && !justSignedOut) {
    const go = new URL("/api/auth/login", url.origin);
    go.searchParams.set("provider", provider);
    go.searchParams.set("next", nextPath);
    return new Response(null, {
      status: 302,
      headers: {
        Location: go.toString(),
        // 60s breaker: cleared by a successful callback.
        "Set-Cookie": `${AUTOTRY_COOKIE}=1; Path=/; HttpOnly; Secure; SameSite=Lax; Max-Age=60`,
        "Cache-Control": "no-store",
      },
    });
  }

  // Already on the sign-in page: serve it rather than redirecting to itself.
  if (onLoginPage) {
    const res = next();
    res.headers.set("Cache-Control", "no-store");
    return res;
  }

  const login = new URL("/login", url.origin);
  login.searchParams.set("next", nextPath);
  return Response.redirect(login, 302);
}
