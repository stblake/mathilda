# Task: Issue #77 — clarify installation (runtime code + libraries), Win/Linux/Mac

**Issue:** Matthias Köppe (SageMath packager): "What do I need to copy into an
installation prefix other than the `Mathilda` executable?"

**Root cause:** the binary is not self-contained (the `src/internal/` `.m` module
tree is loaded at runtime), and there was no `make install` target nor any
documentation of what to ship. The loader already searches an exe-relative FHS
layout (`<exe>/../share/mathilda/internal/`) — nothing populated it.

## Plan

- [x] Add `make install` / `make uninstall` (GNU vars: `DESTDIR`, `PREFIX`/`prefix`,
      `bindir`, `datadir`, `INSTALL*`); installs binary → `$(PREFIX)/bin`, module
      tree → `$(PREFIX)/share/mathilda/internal/`. Lower-case `prefix` default so a
      bare `make` doesn't trip the opt-in `ifdef PREFIX` (`-DMATHILDA_PREFIX`).
- [x] README: new "Installing Mathilda (deploying to a prefix)" section — the two
      things to ship, `make install` + `DESTDIR`/overrides, the 4 resolution
      mechanisms, manual copy, per-platform shared-library story (Win/Linux/Mac).
- [x] README: document the undocumented **PCRE2** (`USE_REGEX`) optional dep
      (prereq list, backends table, apt/brew/dnf install lines).
- [x] Changelog: `## Build & packaging` note in `docs/spec/changelog/2026-09-21.md`
      (no `$VersionNumber` bump — build tooling + prose are contributor-facing).
- [ ] Respond to issue #77 (outward-facing — awaiting go-ahead on commit/push+post).

## Review

**Changes:** `makefile` (+install/uninstall, +.PHONY), `README.md` (install
section + PCRE2), `docs/spec/changelog/2026-09-21.md` (build note). No C source
touched → no version bump/tag.

**Verified (staged install into scratchpad):**
- `make install DESTDIR=<stage> PREFIX=/usr/local` → binary at
  `.../usr/local/bin/Mathilda`, all 11 `.m` files under
  `.../usr/local/share/mathilda/internal/` with `simp/`, `simp/transforms/`,
  `mixed/` preserved.
- Installed binary run from CWD=`/` with `MATHILDA_HOME` unset (proves
  exe-relative resolution + relocatability): `FullSimplify[Sin[x]^2+Cos[x]^2]`→`1`
  (simp tree), `Integrate[1/(x^2+1),x]`→`ArcTan[x]`, `BesselJ[1/2,x]`→elementary
  (bessel.m). No `LoadModule::nofile` on stderr.
- `make uninstall DESTDIR=<stage> PREFIX=/usr/local` → 0 Mathilda files remain.
- Bare `make` does NOT define `-DMATHILDA_PREFIX` (opt-in preserved);
  `make PREFIX=…` does (fallback intact).
