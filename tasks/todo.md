# Task: `Integrate::nonelem` warning for the ParallelMixedTower method

Replicate RischTranscendental's non-elementary warning for the ParallelMixedTower method.
Full plan: `~/.claude/plans/the-method-rischtranscendental-issues-precious-cookie.md`.

## Implementation checklist

- [x] `integrate.c`: add `#include "print.h"` (for `expr_to_string`)
- [x] `integrate.c`: add file-static de-dup flag `g_integrate_nonelem_announced`
- [x] `integrate.c`: add `integrate_announce_nonelementary(Expr*, Expr*)` helper
- [x] `integrate.h`: declare the helper
- [x] `integrate.c`: reset flag after `g_integrate_depth++` (depth==1)
- [x] `integrate.c`: add `pmt_is_nonelementary_certificate` predicate + call helper in `builtin_integrate_pmt`
- [x] `integrate_risch_transcendental.c`: route inline fprintf through the helper
- [x] `integrate.c`: one docstring line for `Integrate\`ParallelMixedTower`
- [x] `tests/test_parallelmixedtower.c`: certificate + stderr-capture tests
- [x] `src/version.h`: bump 0.165 -> 0.166
- [x] `docs/spec/builtins/calculus.md`: note the certificate
- [x] `docs/spec/changelog/2026-09-21.md`: changelog note
- [x] Build, run tests, manual REPL proof, valgrind spot-check
- [ ] Commit + tag v0.166 (awaiting user go-ahead — on `main`, not auto-committing)

## Review

Implemented exactly as planned. The `.m` worker already computed rigorous
`{"not elementary", …}` certificates (the paper's three guarded `NotElementary`
exits); the only gap was C-side, where the dispatcher discarded every list as a
decline. Added a shared `integrate_announce_nonelementary` helper (used by both
RischTranscendental and ParallelMixedTower) with a per-cascade de-dup flag so the
two methods never double-print, and a `pmt_is_nonelementary_certificate` predicate
that distinguishes `{"not elementary", …}` (warn) from `{"failed", …}` (silent).

**Verified:**
- `Integrate[1/(x Log[x+Sqrt[x^2+1]]), x]` and `Integrate[Tan[Sqrt[x^2+1]], x]`
  → `Integrate::nonelem` on stderr, unevaluated (residue certificates).
- `Exp[x^2]` via ParallelMixedTower → `{failed,…}`, correctly silent.
- Elementary flagship `Log[x+Sqrt[x^2+1]]` → integrates, no message.
- Automatic cascade prints the message exactly once (de-dup).
- `$VersionNumber` → 0.166.
- Suites pass: parallelmixedtower_tests, integrate_risch_transcendental_tests,
  integrate_chebychev_tests, dsolve_tests. `make check-c99` clean.
- valgrind: no leak frame references the new code (`integrate.c` absent from all
  leak stacks); the pre-existing `.m`-package/baseline leaks are unrelated.
