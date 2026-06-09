# Task: Move special functions to src/special_functions/

Modules (13): gamma, loggamma, polygamma, pochhammer, eulergamma, zeta,
stieltjesgamma, bernoullib, eulere, polylog, hypergeopfq, fibonacci, lucas
(.c + .h each = 26 files)

- [ ] 1. `git mv` the 26 files into `src/special_functions/`
- [ ] 2. makefile: add `src/special_functions/*.c` wildcard + `-I./src/special_functions`
- [ ] 3. tests/CMakeLists.txt: update 13 source paths + add include_directories
- [ ] 4. Update doc source-path references (special-functions.md, mathematical-constants.md, SPEC.md layout)
- [ ] 5. Build `make -j` cleanly
- [ ] 6. Build & run a relevant test binary
- [ ] 7. Changelog note (docs/spec/changelog/2026-06-08.md)

## Notes
- All modules registered only in core.c; no other source includes their headers.
- Cross-includes (`#include "gamma.h"`) resolve via new `-I` path — no edits to includes needed.
- site/ doc refs are generated (regen later via `make docs`); changelog history left as-is.

## Review
Done. All 13 modules (26 files) moved via `git mv` (tracked as renames, history
preserved). makefile + tests/CMakeLists.txt updated (wildcard, -I path, 13 test
source paths). No source `#include` or `core_init()` edits needed — cross-includes
resolve via the new include path. Stale `src/*.o` from prior builds removed.
Main binary builds clean (`-Wall -Wextra`); REPL evaluates Zeta/Gamma/Fibonacci/
BernoulliB/PolyLog/LucasL correctly; `gamma_tests` passes. Docs updated: SPEC.md
layout, special-functions.md + mathematical-constants.md source paths, and a
changelog note in docs/spec/changelog/2026-06-08.md. site/ refs left for `make docs`
regeneration; historical changelog entries left as-is.
