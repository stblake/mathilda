# Task: Implement StringQ

`StringQ[expr]` → True if expr is a string, False otherwise.
`StringQ[]` (wrong arity) → `StringQ::argx` message, left unevaluated.

## Plan
- [x] sym_names.h: declare `SYM_StringQ`
- [x] sym_names.c: define + intern `SYM_StringQ`
- [x] core.h: declare `builtin_stringq`
- [x] core.c: implement `builtin_stringq` (arg_error on argc != 1); register + ATTR_PROTECTED
- [x] info.c: docstring (no examples)
- [x] docs/spec/builtins/expression-information.md: document StringQ
- [x] docs/spec/changelog/2026-08-03.md: changelog entry
- [x] tests/test_core.c: `test_stringq`, register in main
- [x] Build main + core_tests, run, check no leaks

## Review
Implemented `StringQ[expr]` following the existing `*Q` predicate pattern
(modelled on `NumberQ`/`IntegerQ`/`AtomQ`).

- **Behaviour** (verified against the spec in the task):
  `StringQ["AbC"]`, `StringQ[""]`, `StringQ["123"]` → `True`;
  `StringQ[123]`, `StringQ[1.5]`, `StringQ[x]`, `StringQ[Pi]`,
  `StringQ[{"a","b"}]` → `False`. Argument is evaluated first
  (`StringQ[If[True,"yes",0]]` → `True`). `StringQ[]` emits
  `StringQ::argx: StringQ called with 0 arguments; 1 argument is expected.`
  and leaves the call unevaluated.
- **Definite boolean**, never symbolic (the `*Q` contract) except the
  documented malformed-arity carve-out, which uses the shared
  `builtin_arg_error` helper.
- **Not Listable** — only `Protected`, matching Mathematica.
- **Verification:** `core_tests` "All core tests passed"; `leaks` = 0 leaks /
  0 bytes over a 20,000-iteration all-paths loop; `make check-c99` exit 0;
  docstring registered and viewable via `Information[StringQ]`.
- **Docs:** `docs/spec/builtins/expression-information.md` (worked examples) +
  changelog entry in `docs/spec/changelog/2026-08-03.md`.

Files touched: `src/core.{c,h}`, `src/sym_names.{c,h}`, `src/info.c`,
`tests/test_core.c`, the two docs files.
