# Access control

This site is gated by **Lab Logger SSO**. There is no password and no user list in this repo.

You may open this site if — and only if — you are a member of the **Zernicka-Goetz Lab** in
[Lab Logger](https://www.lab-logger.com). Adding or removing someone there grants or revokes
access here. Nothing else needs to change.

## How it works

1. `/api/auth/login?provider=google|apple` starts a **server-side PKCE** flow against the same
   Supabase project Lab Logger uses. The PKCE verifier and a CSRF `state` are sealed into a
   short-lived signed httpOnly cookie — never put in the URL.
2. `/api/auth/callback` verifies that cookie and the `state`, then exchanges the auth code at
   `/auth/v1/token?grant_type=pkce`.
3. It asks PostgREST whether the user belongs to an allowed lab, **using the user's own access
   token**, so Postgres RLS makes the decision:
   ```
   GET /rest/v1/lab_members?select=lab_id,role&lab_id=in.(<ALLOWED_LAB_IDS>)&limit=1
   apikey: <anon>   Authorization: Bearer <the user's token>
   ```
   A non-member gets **HTTP 200 with an empty array**, not an error — so the check is
   `Array.isArray(rows) && rows.length > 0` (`isMemberResult`). Testing `Array.isArray(rows)`
   alone would admit *everyone*; `tests/session.test.mjs` asserts both branches.
   A **service_role key is never used** — it would bypass RLS entirely.
4. On success the site mints **its own** session cookie: HMAC-SHA256 signed, httpOnly, Secure,
   SameSite=Lax, carrying `{sub, em, nm, rl, adm, exp}` and nothing else. The Supabase
   access/refresh tokens are discarded — this site never acts on your behalf in Lab Logger.

Identity is matched on the **Supabase user id**, never on email: Apple's Hide My Email returns a
`@privaterelay.appleid.com` address that will not match the same person's Google address.

## Enforcement

`middleware.js` (Vercel Edge Middleware) runs **before any file is served**. Its matcher covers
every content path — the `build_*.py` output under `/data`, the project pages, scripts, styles —
excluding only `/api/auth/*` and Vercel internals.

* HTML navigation, unauthenticated → `302` to `/login?next=…`
* anything else, unauthenticated → `401` (so `fetch()` and `curl` get a clean refusal)
* **fails closed**: with `SESSION_SECRET` unset it returns `503` rather than serving data
* `/admin*` and `/api/admin*` return **404** (not 403) for non-admins, so the console's
  existence is never revealed

Caching: every gated response is `Cache-Control: private` with `Vary: Cookie`, so Vercel's shared
CDN can never hand a cached authenticated response to an unauthenticated request. `/data/*` keeps
`max-age=3600` so *browsers* still cache the large files.

## Zero-click sign-in

This site cannot read Lab Logger's Supabase session (different origin; third-party cookies are
dead), so a bounce through Google/Apple is unavoidable — but that bounce is **silent** when a
provider session already exists. The only real friction is our own provider chooser, which is
skipped whenever the provider is known:

* `?provider=google|apple` on the incoming link (Lab Logger can add this via its
  `src/lib/labResources.ts` `passProvider: true` opt-in — **ask for this site to be added there**);
* otherwise the `surf_provider` cookie remembered from the last successful sign-in.

We **never guess** a provider with no signal. Sending an Apple user to Google can silently
authenticate a *different* account that isn't in the lab, so the fallback is a one-time chooser.

Guards, all of which are load-bearing:

* `surf_auto` — a ~60 s breaker so a *failed* silent attempt falls back to the chooser instead of
  looping;
* `/login` is inside the matcher (bookmarks auto-continue) but never redirects to itself;
* sign-out goes to `/login?signedout=1`, and that flag suppresses auto-continue — otherwise
  signing out would be impossible. The remembered provider is deliberately kept so signing back
  in is still one silent bounce.
* `?next=` accepts only same-site paths starting with a single `/` (`safeNext`), or it becomes an
  open redirect.

## Roles

Admin = listed in `ADMIN_USER_IDS`, **or** your Lab Logger `lab_members.role` is `admin` or `pi`.
Admins get `/admin` (usage analytics + project access) and always see every project.

## Per-project access

Lab membership decides whether you get in at all. `project_access` (Neon) optionally narrows
*which projects* a given member may open — a member with **no row sees every project**.

Keyed on the Supabase user id. An admin can restrict someone **before they have ever signed in**
by email: the row is stored as `pending:<email>` and rebound to the real user id on first
sign-in (`bindAccessRow`), so the limit applies from their very first visit.

The lookup **fails open** on a database error — a Neon hiccup must not lock an authenticated lab
member out of the whole site. The gate that matters (lab membership) has already passed.

> The pre-SSO `access` table (keyed on display name) is left untouched as a historical record and
> is no longer read. Drop it whenever you like.

## Environment variables

| Variable | Required | Notes |
|---|---|---|
| `SUPABASE_URL` | yes | Same project as Lab Logger |
| `SUPABASE_ANON_KEY` | yes | Publishable key — safe to expose; RLS enforces |
| `SESSION_SECRET` | yes | `openssl rand -base64 32`. **Rotating it signs everyone out.** |
| `ALLOWED_LAB_IDS` | yes | Comma-separated lab uuids |
| `SESSION_TTL_HOURS` | no | Default 12, clamped 1–720 |
| `ADMIN_USER_IDS` | no | Comma-separated Supabase user ids |
| `SITE_ORIGIN` | no | Overrides the origin used to build `redirect_to` |
| `DATABASE_URL` | no | Neon, for analytics + project access |

⚠️ **Set these before deploying** — the middleware fails closed, so a deploy without
`SESSION_SECRET` returns 503 to everyone.

⚠️ Also add this origin to **Supabase → Authentication → URL Configuration → Redirect URLs** as
`https://<this-site>/**`, or the callback is rejected.

## Tests

```bash
npm test
```

Covers session sign/verify (tamper, wrong secret, expired, missing `exp`, malformed), the
`?next=` open-redirect guard, **both branches of the membership check**, provider inference, and
the middleware enforcement matrix (401 on data, 302 on HTML, 404 on admin, fail-closed 503,
breaker, sign-out).
