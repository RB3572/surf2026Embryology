// GET /api/auth/login?provider=google|apple&next=/some/path
//
// Starts a PKCE login against the SAME Supabase project Lab Logger uses. The PKCE verifier
// is held in a short-lived signed httpOnly cookie (never in the URL), and Supabase returns to
// /api/auth/callback with an auth code. Access tokens are exchanged server-side and never
// reach the browser.

import {
  PKCE_COOKIE, randomToken, sha256b64url, signToken,
  cookieHeader, requiredEnv, originFor, safeNext,
} from "../../lib/session.mjs";

const PROVIDERS = new Set(["google", "apple"]);

export default async function handler(req, res) {
  const missing = requiredEnv();
  if (missing.length) {
    res.statusCode = 500;
    res.setHeader("content-type", "text/plain; charset=utf-8");
    res.setHeader("Cache-Control", "no-store");
    res.end(`Auth is not configured. Set in Vercel: ${missing.join(", ")}.`);
    return;
  }

  const provider = String(req.query.provider || "");
  if (!PROVIDERS.has(provider)) {
    res.statusCode = 400;
    res.setHeader("Cache-Control", "no-store");
    res.end("unknown provider");
    return;
  }

  const secret = process.env.SESSION_SECRET;
  const origin = originFor(req);
  const next = safeNext(req.query.next);

  // PKCE pair for our leg to Supabase + CSRF state, both sealed into a signed cookie so the
  // callback can't be replayed or cross-site forged.
  const verifier = randomToken(32);
  const challenge = await sha256b64url(verifier);
  const state = randomToken(16);
  const pkce = await signToken(
    { v: verifier, s: state, n: next, p: provider,
      exp: Math.floor(Date.now() / 1000) + 600 },
    secret);

  const authUrl = new URL(`${process.env.SUPABASE_URL}/auth/v1/authorize`);
  authUrl.searchParams.set("provider", provider);
  authUrl.searchParams.set("redirect_to", `${origin}/api/auth/callback?state=${state}`);
  authUrl.searchParams.set("code_challenge", challenge);
  authUrl.searchParams.set("code_challenge_method", "s256");

  res.setHeader("Set-Cookie", cookieHeader(PKCE_COOKIE, pkce, 600));
  res.setHeader("Cache-Control", "no-store");
  res.statusCode = 302;
  res.setHeader("Location", authUrl.toString());
  res.end();
}
