// GET|POST /api/auth/logout — drop this site's session cookie.
// This signs the user out of THIS SITE only; their Lab Logger / Google / Apple session is
// untouched, so signing back in is one click.

import {
  SESSION_COOKIE, PKCE_COOKIE, AUTOTRY_COOKIE, clearCookieHeader,
} from "../../lib/session.mjs";

export default function handler(req, res) {
  // NOTE: surf_provider is deliberately KEPT, so signing back in later is a single silent
  // bounce rather than another chooser. The ?signedout=1 guard in middleware is what stops us
  // re-authenticating them immediately — without it, signing out would be impossible.
  res.setHeader("Set-Cookie", [
    clearCookieHeader(SESSION_COOKIE),
    clearCookieHeader(PKCE_COOKIE),
    clearCookieHeader(AUTOTRY_COOKIE),
  ]);
  res.setHeader("Cache-Control", "no-store");
  res.statusCode = 302;
  res.setHeader("Location", "/login?signedout=1");
  res.end();
}
