# Task: Implement PadRight in list.c

PadRight is the right-padding mirror of WL PadLeft. Implement a comprehensive,
multidimensional engine.

## Plan
- [x] sym_names.h / sym_names.c: register `SYM_PadRight`
- [x] list.h: declare `builtin_padright`
- [x] list.c: implement `builtin_padright` + helpers
      - [x] `pr_pad_at` — cyclic index into padding block by coord path
      - [x] `pr_full_dims` — max dims of a ragged array (for `PadRight[list]` / Automatic)
      - [x] `pr_build` — recursive multidim padding engine
      - [x] argb validation (1..4 args), head preservation
- [x] list_init: register builtin + ATTR_PROTECTED
- [x] info.c: docstring
- [x] tests/test_padright.c + CMakeLists entry
- [ ] docs/spec/builtins/<list>.md + changelog
- [ ] build, run tests, valgrind

## Semantics (PadRight = mirror of PadLeft)
- PadRight[list, n]            pad with 0 on the right to length n
- PadRight[list, n, x]         pad with element x
- PadRight[list, n, {x...}]    cyclic padding (first list elem -> P[0])
- PadRight[list, n, pad, m]    margin m of padding on the LEFT
- PadRight[list, n, pad, -m]   truncates first m elements
- negative n                   pads on the LEFT instead
- PadRight[list, {n1,n2,...}]  nested full array, dim n_i at level i
- PadRight[list]               pad ragged array to full with 0
- PadRight[list, Automatic, x] pad to full array with x
- head need not be List (preserved)

## Review
Done. Implemented BOTH PadRight and PadLeft (user follow-up) in src/list.c as
exact mirrors sharing one recursive engine `pr_build(... pad_left ...)`.

- Engine: pr_pad_at (cyclic block index), pr_scan_dims (full dims), pr_build
  (recursive multidim layout), pad_dispatch (arg parsing, shared by both).
- Key correctness point: padding phase is data-anchored. For PadLeft the
  `list_start + L` anchor makes L cancel (= N - m), so empty pure-padding
  sub-rows in a margined multidim pad share the data rows' tiling. Only a
  wholly-empty *top-level* list anchors at position 0 (empty_anchor flag).
  This made the documented {5,5},{{x,y},{z}},{1,2} example match WL exactly.
- All documented WL examples reproduced exactly (verified in REPL).
- Symbols in sym_names.{c,h}; docstrings in info.c; ATTR_PROTECTED in list_init.
- Tests: tests/test_padright.c + tests/test_padleft.c (19 cases each), both
  registered in tests/CMakeLists.txt. All pass; list_tests still pass.
- Valgrind: only the ~12,800B/400-block macOS dyld baseline; no list.c frames
  in any leak stack -> no new leaks.
- Docs: docs/spec/builtins/structural-manipulation.md + changelog 2026-06-15.md.
