# Remove the `Locked` attribute (`ATTR_LOCKED`)

Full removal — nothing in open-source Mathilda should be permanently locked.
`Remove`'s self-protection survives via `Protected`.

## A. C core — flag + machinery
- [ ] `src/attr.h` — delete `ATTR_LOCKED` define, leave gap comment
- [ ] `src/attr.c` — set_attributes guard, set/clear_attributes_for_symbol guards, string_to_attribute mapping, builtin_attributes count+emit, scan attr_init
- [ ] `src/core.c` — Remove seeding; Unset/ClearAll/Remove masks; Protect/Unprotect early-returns
- [ ] `src/eval.c` — value-list assign, UpSet, TagSet masks
- [ ] `src/options_builtin.c` — delete SetOptions::locked guard
- [ ] `src/options.h` — fix comment
- [ ] `src/sym_names.{c,h}` — remove SYM_Locked decl/def/intern

## B. C docstrings
- [ ] `src/info.c` — SetOptions, ClearAll, Remove, Protect, Unprotect

## C. Authored docs
- [ ] `SPEC.md` attribute table
- [ ] `docs/spec/builtins/assignment-and-rules.md`
- [ ] `docs/spec/builtins/expression-information.md`
- [ ] `docs/spec/changelog/2026-10-05.md` — add removal entry + reconcile same-week
- [ ] `book/chapters/10-internals.tex` caption

## D. Tests
- [ ] `tests/test_clearall_remove_protect.c` — Attributes[Remove] literal
- [ ] `tests/test_upvalues.c` — convert Locked case to Protected

## E. Generated doc pages (hand-edit)
- [ ] `frontend/public/refpages/` assignment-and-rules + expression-information
- [ ] `site/impl/`, `site/overlays/`, `site/docs/documentation/`

## F. Version + tag
- [ ] `src/version.h` bump to 0.301 (tag on commit, when asked)

## Verification
- [ ] Build clean; check-messages / check-c99
- [ ] Final greps zero
- [ ] Unit tests (incl. -DUP_VALUES)
- [ ] REPL behavior checks
- [ ] Rebuild code-review graph

## Review

**Done — `ATTR_LOCKED` fully removed (v0.301).**

C core: flag deleted from `attr.h` (gap comment left, matching the bit-11
precedent); `SYM_Locked` removed from `sym_names.{c,h}`; the `"Locked"` name
mapping + `Attributes[]` count/emit removed from `attr.c`; the three attribute-
setter guards dropped; `Remove` no longer seeded Locked (`core.c`); all combined
`(ATTR_PROTECTED | ATTR_LOCKED)` guards across `core.c`/`eval.c` reduced to
`ATTR_PROTECTED`; `Protect`/`Unprotect` Locked early-returns deleted; the
`SetOptions::locked` guard + message removed (`options_builtin.c`); comments in
`options.h` fixed. Docstrings updated in `info.c`.

Docs: `SPEC.md` attribute table, `docs/spec/builtins/{assignment-and-rules,
expression-information}.md`, `book/chapters/10-internals.tex` caption, and the
current-week changelog (`2026-10-05.md`) with a removal entry + same-week
reconciliation. 19 generated pages under `frontend/public/refpages/` and `site/`
hand-edited. Historical (pre-this-week) changelogs left as dated records.

Tests: `tests/test_clearall_remove_protect.c` (`Attributes[Remove]` →
`{HoldAll, Protected}`) and `tests/test_upvalues.c` (`test_upset_locked` →
`test_upset_protected`). Both suites green.

Verification: clean `make` build; `make check-messages` and `make check-c99`
pass; `grep ATTR_LOCKED\|SYM_Locked src/` → only the intentional gap comment;
repo-wide `\bLocked\b` audit → only intended + coincidental hits. REPL:
`Attributes[Remove]` = `{HoldAll, Protected}`; `SetAttributes[x, Locked]` is an
ignored no-op (`{}`); a user symbol is Protect→Unprotect→Remove-able;
`$VersionNumber` = 0.301.

Not committed/tagged yet (awaiting user request). On commit: tag `v0.301`.
