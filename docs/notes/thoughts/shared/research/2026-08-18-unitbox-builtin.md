---
date: 2026-08-18T14:15:54-07:00
researcher: Michael Sollami
git_commit: dfaa8aa1df0518dbf7613151196a476411bd2baa
branch: feature/unitbox
repository: mathilda
topic: "Adding UnitBox as a new built-in (Ramp/UnitStep as siblings)"
tags: [research, codebase, piecewise, unitbox, ramp, unitstep, attributes, docs-spec, tests]
status: complete
last_updated: 2026-08-18
last_updated_by: Michael Sollami
last_updated_note: "Added Decisions section resolving the four Open Questions"
---

# Research: Adding UnitBox as a new built-in

**Date**: 2026-08-18T14:15:54-07:00
**Researcher**: Michael Sollami
**Git Commit**: dfaa8aa1df0518dbf7613151196a476411bd2baa
**Branch**: feature/unitbox
**Repository**: mathilda

## Research Question

How are `Ramp` and `UnitStep` declared, registered, and implemented, can `UnitBox`
reuse `ustep_class`, what attributes does each set, how are the `docs/spec/builtins`
entries structured, what does `tests/test_clip.c` look like structurally, and does
`UnitBox` already exist anywhere in the tree?

## Summary

`Ramp` and `UnitStep` are both implemented in `src/piecewise.c`, registered in
`piecewise_init()`, documented via `symtab_set_docstring()` in `src/info.c`, and
spec'd in `docs/spec/builtins/elementary-functions.md`. Both lean on a shared
sign-classification helper, `ustep_class()` (`src/piecewise.c:491`), which resolves
an argument's sign against zero (exact, symbolic-certified, or unknown). `UnitBox`
needs a **band** test (`-1/2 <= x <= 1/2`), which `ustep_class` alone cannot answer
since it only classifies sign against a single threshold of 0 — `UnitBox[x]` would
need two such comparisons (`x >= -1/2` and `x <= 1/2`), most naturally built by
calling `ustep_class` twice on shifted arguments (`x + 1/2` and `1/2 - x`), or by
writing an analogous three-way classifier. `UnitBox` does not exist anywhere in the
repository (confirmed by exhaustive grep). `tests/test_clip.c` is the closest test
exemplar: individual `static void test_*(void)` functions driven through
`assert_eval_eq`/`assert_eval_startswith` (string-compare on printed output), run
from `main()` via the `TEST(...)` macro, plus a dedicated memory-loop test for
valgrind coverage. Its CMake target (`clip_tests`) is built but **not** wired into
`ctest` via `add_test()` — worth checking whether a new `unitbox_tests` target
should follow the same (registered or unregistered) pattern as its neighbors.

## Detailed Findings

### Registration (`piecewise_init`)

- `src/piecewise.c:22-45` — `piecewise_init()` registers all of `Floor`, `Ceiling`,
  `Round`, `IntegerPart`, `FractionalPart`, `UnitStep`, `Ramp` via
  `symtab_add_builtin`, then sets attributes in a second pass.
- `src/piecewise.c:28` — `symtab_add_builtin("UnitStep", builtin_unitstep);`
- `src/piecewise.c:29` — `symtab_add_builtin("Ramp", builtin_ramp);`

### Attributes — in code

- `src/piecewise.c:36-39`:
  ```c
  /* UnitStep: Listable threads over lists, NumericFunction marks it for N,
   * Orderless because UnitStep[x1,...,xn] is symmetric in its arguments. */
  symtab_get_def("UnitStep")->attributes |=
      (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE | ATTR_ORDERLESS);
  ```
  → `UnitStep`: `Protected | NumericFunction | Listable | Orderless`.
- `src/piecewise.c:41-44`:
  ```c
  /* Ramp: Listable, NumericFunction, Protected -- Mathematica's attribute
   * set exactly. NOT Orderless: Ramp is unary. */
  symtab_get_def("Ramp")->attributes |=
      (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
  ```
  → `Ramp`: `Protected | NumericFunction | Listable` (explicitly **not** Orderless —
  called out in the comment because it's unary).
- `UnitBox[x]` is unary like `Ramp`, so by the same logic it should get
  `Protected | NumericFunction | Listable` and **not** `Orderless`.

### Attributes — in the spec doc

- `docs/spec/builtins/elementary-functions.md:584-585`:
  ```
  **Features**:
  - `Listable`, `NumericFunction`, `Orderless`, `Protected`.
  ```
  (UnitStep entry — matches the code exactly, just reordered alphabetically-ish.)
- `docs/spec/builtins/elementary-functions.md:629-630`:
  ```
  **Features**:
  - `Listable`, `NumericFunction`, `Protected`.
  ```
  (Ramp entry — matches the code; no `Orderless`.)

### `ustep_class` — what it does and reuse for UnitBox

- `src/piecewise.c:407` — `enum { USTEP_UNKNOWN = -1, USTEP_NEG = 0, USTEP_NONNEG = 1 };`
- `src/piecewise.c:413-443` — `ustep_real_sign(Expr* a)`: returns -1/0/+1 for a
  concrete real numeric leaf (`EXPR_INTEGER`, `EXPR_BIGINT`, `EXPR_REAL`,
  `EXPR_MPFR` under `USE_MPFR`), or exact int64/GMP rationals, or the sentinel
  `-2` if `a` isn't a pure real number.
- `src/piecewise.c:447-449` — `ustep_is_pos_infinity`: true only for the bare
  `Infinity` symbol (not `ComplexInfinity`/`DirectedInfinity`).
- `src/piecewise.c:457-488` — `ustep_certify(Expr* x)`: for an exact symbolic real
  (`Pi`, `Sqrt[2]`, ...) that no leaf branch resolved, numericalizes at increasing
  MPFR precision (`64` up to `1<<14` bits, `*4` per step,
  `src/piecewise.c:459-461`) and accepts a sign only once two successive
  precisions agree (mirrors `$MaxExtraPrecision` — never guesses). Falls back to a
  single machine-precision probe when `USE_MPFR` is undefined
  (`src/piecewise.c:479-487`).
- `src/piecewise.c:491-511` — `ustep_class(Expr* x)`: the public classifier.
  1. `+Infinity` → `USTEP_NONNEG` (`:492`).
  2. `-Infinity` (via `is_minus_infinity`, `src/piecewise.c:90-99`) → `USTEP_NEG`
     (`:493`).
  3. `ustep_real_sign` on `x` directly; if resolved (`!= -2`), returns
     `NEG`/`NONNEG` treating `0` as `NONNEG` (`:495-496`).
  4. If `x` is `Complex[re, im]` (via `is_complex`), a non-zero imaginary part →
     `USTEP_UNKNOWN`; a zero imaginary part falls through to classify the real
     part (`:498-507`).
  5. Otherwise, `ustep_certify(x)` for exact symbolic reals (`:510`).

**Can `UnitBox` reuse it?** Only as a building block, not directly — `ustep_class`
answers a single question (sign of `x` against `0`). `UnitBox[x] = 1` iff
`-1/2 <= x <= 1/2`, i.e. two threshold tests. The natural reuse is to call
`ustep_class` twice on shifted arguments:
- `ustep_class(x + 1/2)` — `NEG` means `x < -1/2` → `UnitBox` is `0`.
- `ustep_class(1/2 - x)` — `NEG` means `x > 1/2` → `UnitBox` is `0`.
- If neither call is `NEG` (both `NONNEG`), `UnitBox` is `1` (mirroring `UnitStep`'s
  convention that the value *at* the boundary, `x = ±1/2`, is `1`, matching
  Mathematica's closed interval `[-1/2, 1/2]`).
- If either call returns `USTEP_UNKNOWN`, `UnitBox` is left unevaluated (`NULL`) —
  same policy as `Ramp`/`UnitStep` for undecidable signs.

This reuses all of `ustep_class`'s exact/rational/MPFR/certification machinery
for free, at the cost of building two shifted subexpressions (`Plus[x, 1/2]`,
`Plus[Rational[1,2], Times[-1,x]]` or similar) per call — cheap, and avoids
duplicating `ustep_real_sign`/`ustep_certify`.

### How `Ramp` preserves argument precision/numeric head

- `src/piecewise.c:562-577` (comment) — explains the design rationale: the zero
  returned for a negative argument carries the *argument's own* exactness, so a
  buffer of one numeric head maps to a result of the same head (`Ramp[-1.]` is
  `0.`, `Ramp[-3]` is exact `0`) — this is what makes a uniform buffer answer
  legal, unlike `Clip` which mixes heads when a Real value clips to an exact bound.
- `src/piecewise.c:578-593` — `builtin_ramp`:
  ```c
  Expr* builtin_ramp(Expr* res) {
      if (res->type != EXPR_FUNCTION) return NULL;
      if (res->data.function.arg_count != 1) return NULL;

      Expr* x = res->data.function.args[0];
      int cls = ustep_class(x);
      if (cls == USTEP_UNKNOWN) return NULL;
      if (cls == USTEP_NONNEG) return expr_copy(x);

      /* Negative: zero, at the argument's own precision. */
      if (x->type == EXPR_REAL) return expr_new_real(0.0);
  #ifdef USE_MPFR
      if (x->type == EXPR_MPFR) return expr_new_mpfr_bits(mpfr_get_prec(x->data.mpfr));
  #endif
      return expr_new_integer(0);
  }
  ```
  - Non-negative branch: returns `expr_copy(x)` verbatim (`:585`) — the argument's
    exact type/value/precision is untouched.
  - Negative branch (`:588-592`): dispatches on `x->type` to build a same-headed
    zero: `EXPR_REAL` → `expr_new_real(0.0)`; `EXPR_MPFR` (guarded by `USE_MPFR`) →
    `expr_new_mpfr_bits(mpfr_get_prec(x->data.mpfr))` (a zero at the *same bit
    precision* as `x`); otherwise → exact `expr_new_integer(0)`.
  - For `UnitBox`, the return values are always the *constants* `0` or `1`
    (never the argument itself), so this precision-preservation pattern from
    `Ramp` does not directly transfer — `UnitBox` is closer to `UnitStep` in that
    respect (`UnitStep` always returns exact `0`/`1`, per
    `src/piecewise.c:521` comment and `docs/spec/builtins/elementary-functions.md:586-588`).

### Docstrings (`symtab_set_docstring`)

- `src/info.c:2791-2799`:
  ```c
  symtab_set_docstring("UnitStep",
      "UnitStep[x]\n"
      "\tgives 0 for x < 0 and 1 for x >= 0 (the value at 0 is 1).\n"
      "UnitStep[x1, x2, ...]\n"
      "\tgives 1 only when none of the xi are negative, otherwise 0.\n"
      "UnitStep[] is 1. The result is always exact. Exact symbolic real\n"
      "arguments are resolved by numerical certification; non-real or\n"
      "unresolved arguments are left unevaluated. UnitStep is Listable and\n"
      "Orderless.");
  ```
- `src/info.c:2800-2809`:
  ```c
  symtab_set_docstring("Ramp",
      "Ramp[x]\n"
      "\tgives x for x >= 0 and 0 for x < 0 -- the positive part of x, and\n"
      "\tthe standard spelling of a rectified linear unit.\n"
      "\n"
      "The zero returned for a negative argument carries the argument's own\n"
      "exactness: Ramp[-1.] is 0. and Ramp[-3] is the exact 0, so a Real\n"
      "list maps to a Real list and an integer list to an integer one.\n"
      "Non-real arguments, and symbolic ones whose sign cannot be decided,\n"
      "are left unevaluated. Ramp is Listable and a NumericFunction.");
  ```
- Both are registered in `src/info.c`, not `src/piecewise.c` — docstrings live in a
  central `info.c` module regardless of where the builtin itself is implemented.

### Interning (`sym_names.c`)

- `src/sym_names.c:1429` — `SYM_UnitStep = intern_symbol("UnitStep");` — `UnitStep`
  has an interned symbol constant used for fast identity comparisons elsewhere
  (e.g. `src/piecewise.c:557`, `src/expand_power.c:226`,
  `src/calculus/deriv.c:1838`). `Ramp` has **no** corresponding `SYM_Ramp` entry —
  grep of `src/` shows `Ramp` referenced only by its string literal
  (`"Ramp"`, `src/piecewise.c:29,43`), never as an interned symbol constant. A new
  `UnitBox` would need a `SYM_UnitBox` only if some other module needs to test the
  head by identity (as `UnitStep` is, for its derivative/series/compile special
  cases); otherwise the `Ramp` precedent (string-literal registration only) is the
  simpler baseline.

### Other call sites that know about `UnitStep` (not `Ramp`)

`UnitStep` is treated specially in several places `Ramp` is not — worth checking
whether `UnitBox` needs the same treatment eventually (out of scope for a first
pass, but flagged for completeness):
- `src/pack.c:1016` — `UnitStep` (not `Ramp`) is in the packed-array `AWARE` list
  (int64 buffer fast path): `"Floor", "Ceiling", "Round", "IntegerPart", "Sign", "UnitStep",`.
  Comment at `src/pack.c:1011-1015` explains why: these are all "trivial
  comparisons" on an int64 buffer. `Ramp` is conspicuously absent from this list
  (confirmed: `grep -n "\"Ramp\"" src/pack.c` → no output) — it has its own buffer
  path documented separately in `docs/spec/builtins/elementary-functions.md:641-643`
  ("packed-arrays.md"), not through this int64 `AWARE` mechanism.
- `src/compile/compile_infer.c:309` and `src/compile/compile_emit_mathfn.c:48,202`
  — `Compile[]` type inference/codegen special-cases `UnitStep` by name.
- `src/calculus/limit.c:265`, `src/calculus/series.c:2412,2951`,
  `src/calculus/integrate_newton_leibniz.c:128`, `src/calculus/deriv.c:1838` —
  calculus modules (`Limit`, `Series`, `Integrate`, `D`) special-case `UnitStep` by
  name/symbol, e.g. the derivative rule at `deriv.c:1838` builds
  `Piecewise[{{Indeterminate, xi == 0}}, 0]` factors (matches the spec doc's
  "Derivative" section, `docs/spec/builtins/elementary-functions.md:600-601,616-620`).
- None of these mention `Ramp`, `Clip`, or any other piecewise sibling — this
  special-casing is UnitStep-specific, presumably because it's the more
  fundamental Heaviside primitive that derivative/series/compile logic composes
  with. A first `UnitBox` implementation does not need to touch any of these
  files; it's pure scope-creep for the initial builtin.

### `docs/spec/builtins/elementary-functions.md` — entry structure & insertion point

Verbatim entries:

- **`UnitStep`** — `docs/spec/builtins/elementary-functions.md:575-621` (heading,
  two description paragraphs, `**Features**:` bullet list, `**Derivative**` note,
  closing ` ```mathematica ``` ` example block with `In[n]:=`/`Out[n]=` pairs).
- **`Ramp`** — `docs/spec/builtins/elementary-functions.md:623-660` (same shape:
  heading, description, `**Features**:` bullets, closing example block).

Template (also cross-checked against `Chop`/`Clip` neighbors):
1. `## <FunctionName>` — level-2 heading, exact head name.
2. Prose description paragraph(s), call forms in backticks, LaTeX-style `$...$`
   inline math for inequalities.
3. `**Features**:` bolded label + `-` bullet list: attributes first
   (`` `Listable`, `NumericFunction`, ... ``), then behavior notes as sub-bullets.
4. Optional bolded inline label for extra sections (e.g. `**Derivative**` for
   `UnitStep`; `Clip`/`Chop` have their own optional labels like
   `**Complex handling**`).
5. Closing ` ```mathematica ... ``` ` fenced block of `In[n]:=`/`Out[n]=` pairs.
6. One blank line, then the next `## ` heading.

**Insertion point**: immediately **before** `## UnitStep` at
`docs/spec/builtins/elementary-functions.md:575` (after the blank line at `:574`,
which follows the prior "Piecewise and Rounding Functions" entry's closing fence
at `:573`). This is consistent both alphabetically (`UnitBox` < `UnitStep`) and
thematically (grouped with the other Heaviside/ramp-family entries).

### Changelog

- Latest weekly file: `docs/spec/changelog/2026-08-17.md` (Monday of the current
  ISO week, per `CLAUDE.md`/`SPEC.md` convention).
- Heading format — `docs/spec/changelog/2026-08-17.md:1-3`:
  ```
  # Changelog: week of 2026-08-17 (Mon) – 2026-08-23 (Sun)

  Feature additions and fixes recorded during this week.
  ```
- Example builtin-addition entry shape —
  `docs/spec/changelog/2026-08-17.md:909-915` (`## Vector analysis: Grad, Div,
  Curl, Laplacian (2026-08-17)` heading, prose paragraph naming the new module/
  file, attribute/registration notes, then `-` bullets).
- A `UnitBox` changelog entry would follow the same pattern:
  `## UnitBox (2026-08-18)` heading (or similar), noting the file
  (`src/piecewise.c`), registration, attributes, and docstring location.

### `Mathilda_spec.md` index

- `Mathilda_spec.md:23` — `Each category lives in [`docs/spec/builtins/`](docs/spec/builtins/):`
- `Mathilda_spec.md:25` — table header `| Category | File |`.
- `Mathilda_spec.md:43` — the relevant category row:
  ```
  | Elementary functions (trig, hyperbolic, log, exp) | [`builtins/elementary-functions.md`](docs/spec/builtins/elementary-functions.md) |
  ```
  This row's description text doesn't enumerate individual functions (unlike a
  few other rows), so adding `UnitBox` requires **no edit** to `Mathilda_spec.md`
  itself — only the changelog table entry described in
  `Mathilda_spec.md:89`: "New entries land in the file for the current week... a
  change touches a built-in's documented behavior, the corresponding
  `docs/spec/builtins/*.md` file is updated as well."

### `tests/test_clip.c` — structural exemplar

- Includes: `tests/test_clip.c:16-28` — `core.h`, `eval.h`, `expr.h`, `parse.h`,
  `print.h`, `symtab.h`, `test_utils.h`, plus `<gmp.h>`, `<stdbool.h>`, `<stdio.h>`,
  `<stdlib.h>`, `<string.h>`.
- Harness macros come from `tests/test_utils.h`, not libc `<assert.h>`:
  - `assert_eval_eq()` — `tests/test_utils.h:18-30` (parses, evaluates, compares
    printed string via `strcmp`, frees both trees).
  - `assert_eval_startswith()` — `tests/test_utils.h:35-50` (`strncmp` prefix
    check, same free pattern).
  - `TEST(name)` macro — `tests/test_utils.h:52` — prints "Running test: %s" then
    calls `name()`.
  - `ASSERT`/`ASSERT_MSG` — `tests/test_utils.h:74-80`, `:93-99` — `exit(1)` on
    failure, deliberately not the libc `assert()` (comment at
    `tests/test_utils.h:54-57` explains: CMake `Release` builds pass `-DNDEBUG`,
    which would silently no-op a bare `assert()`).
- Structure: one `static void test_*(void)` per behavior/case group — 30 such
  functions in `tests/test_clip.c` (lines 34, 40, 45, 57, 65, 71, 82, 91, 97, 106,
  111, 115, 120, 131, 138, 148, 153, 158, 163, 174, 179, 189, 198, 207, 219, 223,
  228, 235, 244, 261, 297), covering: default interval, explicit
  `{min,max}`/replacement bounds, symbolic constants (Pi/E), ±Infinity, list
  threading (incl. nested), complex rejection, symbolic pass-through, arbitrary
  precision (N / MPFR), argument-shape edge cases, attributes, a dedicated
  memory-hygiene loop, and a packed/NDArray three-way scan.
- `main()` — `tests/test_clip.c:317-379` — calls `symtab_init(); core_init();`
  (`:318-319`), then every test via `TEST(...)` in a fixed, comment-banner-grouped
  order, prints `"All clip_tests passed.\n"` (`:377`), `return 0;`.
- Expressions are built via `parse_expression()` and evaluated via `evaluate()`
  inside the macros; comparison is by **printed-string equality** (`expr_to_string`
  + `strcmp`/`strncmp`), not `expr_eq()` or manual `expr_new_*` construction. Two
  test functions (`test_clip_attributes` at `:244`, `test_clip_memory_loop` at
  `:261`) call `parse_expression`/`evaluate`/`expr_to_string` directly instead of
  through the macros.
- Memory pattern: `assert_eval_eq`/`assert_eval_startswith` free both the parsed
  and evaluated tree internally (`test_utils.h:22,29,39,49`);
  `test_clip_memory_loop` (`tests/test_clip.c:261-291`) runs 20 reps over 14 fixed
  Clip expressions, freeing both trees every inner iteration specifically "so
  valgrind has many chances to catch a leak in any path" (comment,
  `tests/test_clip.c:262-264`).
- CMake registration — `tests/CMakeLists.txt:1276-1277`:
  ```
  add_executable(clip_tests test_clip.c $<TARGET_OBJECTS:mathilda_common>)
  target_include_directories(clip_tests PRIVATE ../src)
  ```
  **No** `add_test(NAME clip_tests COMMAND clip_tests)` line exists anywhere in
  `tests/CMakeLists.txt` for this target — it builds but is invisible to `ctest`
  (per the file's own top-of-file warning, `tests/CMakeLists.txt:4-7`), unlike
  several sibling targets that do have a matching `add_test(...)` (e.g.
  `core_algebra_tests` at `:826`, `graph_tests` at `:849`). `make check-tests`
  (`tools/run_test_suite.sh`) is the mechanism described in `SPEC.md` §9 that runs
  *every* built binary regardless of `ctest` registration, so a `UnitBox` test
  binary would still be exercised by `make check-tests` even if left off
  `add_test()` like `clip_tests` is — but adding the `add_test()` line would be
  the more complete choice for a new test file.

### Confirmed: `UnitBox` does not exist anywhere in the repo

Exhaustive `grep -rn`/`grep -rni` for `"UnitBox"`/`"unitbox"` across the whole
repository (`src/`, `tests/`, `docs/`, `docs/spec/`, `Mathilda_spec.md`,
`src/internal/*.m`, `makefile`, `makefile_tests`, `tests/CMakeLists.txt`, and every
hidden dir) returns **zero content matches**. The only occurrences of the string
"unitbox" anywhere are git metadata for the current branch itself
(`.git/refs/heads/feature/unitbox`, `.git/logs/refs/heads/feature/unitbox`) — not
file content.

## Code References

- `src/piecewise.c:22-45` — `piecewise_init()`: registration + attributes for all
  seven piecewise builtins.
- `src/piecewise.c:407-511` — `ustep_class` and its helpers
  (`ustep_real_sign`, `ustep_is_pos_infinity`, `ustep_certify`).
- `src/piecewise.c:523-560` — `builtin_unitstep`.
- `src/piecewise.c:562-593` — `builtin_ramp` (with the precision-preservation
  design comment).
- `src/info.c:2791-2809` — docstrings for `UnitStep` and `Ramp`.
- `src/sym_names.c:1429` — `SYM_UnitStep` interned symbol (no `SYM_Ramp`
  equivalent exists).
- `src/pack.c:1011-1029` — packed-array `AWARE` list; `UnitStep` is in it, `Ramp`
  is not.
- `docs/spec/builtins/elementary-functions.md:575-660` — `UnitStep`/`Ramp` spec
  entries, and the insertion point for `UnitBox` immediately before `:575`.
- `docs/spec/changelog/2026-08-17.md:1-3,909-915` — current week's changelog file
  heading and an example builtin-addition entry.
- `Mathilda_spec.md:23,25,43,89` — category index row and changelog-update
  convention note.
- `tests/test_clip.c` (whole file) and `tests/test_utils.h:18-99` — closest test
  exemplar and its shared harness macros.
- `tests/CMakeLists.txt:1276-1277` — `clip_tests` target registration (unregistered
  with `ctest`).

## Architecture Insights

- **Shared classification helpers over duplicated sign logic.** `UnitStep` and
  `Ramp` both reduce to the same `ustep_class` primitive; a `UnitBox`
  implementation should follow this pattern by composing `ustep_class` on shifted
  arguments (`x + 1/2`, `1/2 - x`) rather than re-deriving sign/certification
  logic from scratch — this also gives it exact-rational and MPFR-certified
  symbolic handling for free.
- **Return-value exactness policy varies by function shape.** `UnitStep` and
  (implicitly) `UnitBox` return only the constants `0`/`1`, always exact,
  regardless of input precision. `Ramp` is different: it returns the argument
  itself unchanged in the non-negative case, so its "exactness preservation" is a
  non-issue for `UnitBox`, which never returns its argument.
- **Docstrings and attributes are set in different files than the implementation.**
  Implementation lives in `src/piecewise.c`; docstrings live centrally in
  `src/info.c`; symbol interning (`sym_names.c`) is opt-in and only needed when
  another module must test the head by identity (as several calculus/compile
  modules do for `UnitStep`, but none do for `Ramp`).
- **The packed-array `AWARE` list and calculus/compile special-casing are
  independent opt-ins**, not something every new piecewise builtin automatically
  gets — `Ramp` shipped without either, and it's a legitimate, complete builtin.
  A first `UnitBox` pass can reasonably follow `Ramp`'s minimal footprint rather
  than `UnitStep`'s more deeply wired one.
- **Test registration is inconsistent by design tolerance, not oversight.**
  `clip_tests` builds but isn't in `add_test()`; the project's own `make
  check-tests` tooling exists precisely because `ctest` registration is not a
  reliable coverage signal in this repo.

## Historical Context (from thoughts/)

None found — this repository had no pre-existing `thoughts/` directory prior to
this research document (created it at `thoughts/shared/research/` to hold this
file, per the research-command convention).

## Related Research

None found.

## Open Questions

- Should `UnitBox` get a `SYM_UnitBox` interned constant in `sym_names.c`? Not
  needed for a `Ramp`-shaped minimal implementation, but would be needed if
  `UnitBox` is later special-cased in `D`/`Series`/`Limit`/`Compile[]` the way
  `UnitStep` is.
- Should `UnitBox` be added to the `src/pack.c:1016` `AWARE` list for the packed
  int64 buffer fast path? `UnitStep` is there; `Ramp` isn't; `UnitBox`'s
  boolean-valued, comparison-only nature makes it look more like `UnitStep` than
  `Ramp` in this respect, but that's a separate, addable-later concern.
- Exact boundary convention: this research assumes `UnitBox[±1/2] = 1` (closed
  interval, matching Mathematica and the analogous `UnitStep[0] = 1` convention)
  — worth confirming against Mathematica's actual `UnitBox` docs before
  implementing, since the user's own definition in the prompt ("`1` for
  `-1/2 <= x <= 1/2`") already confirms this.
- Whether a new `tests/test_unitbox.c` target should include the `add_test()`
  line that `clip_tests` lacks (recommended: yes, for stricter CI coverage than
  the existing precedent).

## Decisions

The four open questions above were resolved by explicit user ruling on
2026-08-18. Recorded here as decisions, with reasoning, so a session with no
prior history can act on them directly.

1. **No `SYM_UnitBox` interned symbol constant.** There are 871 interned
   constants in `sym_names.c`, and they exist specifically for symbols that
   need to be compared by pointer identity in a hot path — e.g. `SYM_UnitStep`
   exists because `D`, `Series`, and `Compile[]` special-case `UnitStep` by
   name. `Ramp` has no interned constant and works fine as a plain builtin.
   Adding one for `UnitBox` with nothing to compare it against would be dead
   weight. This was conditional on one check, not an assumption: if the
   implementation plan wants `UnitBox` to flow through
   `piecewise_interval()`/`interval_apply_function()` the way `Floor`/`Ceiling`
   do, a symbol would be needed for that path. That check has been done —
   `interval_apply_function()` (`interval.c:906-995`) dispatches by `strcmp`
   on the head name, not by symbol pointer, so even routing through it would
   not itself require an interned constant. Separately, neither `UnitStep` nor
   `Ramp` — `UnitBox`'s two direct siblings — currently thread through
   `piecewise_interval()` at all; only `Floor`/`Ceiling` do, because they are
   monotone, and `interval_apply_function`'s dispatch table has no entry for
   `UnitStep` or `Ramp` either. `UnitBox` is a non-monotone, two-sided box
   function, so it does not fit that monotone/region model without new
   interval-arithmetic code. Decision: `UnitBox` will NOT thread through
   `piecewise_interval()` in this first implementation, consistent with its
   siblings — and therefore no `SYM_UnitBox` is needed, for this or any other
   reason currently in scope.
2. **No `src/pack.c` `AWARE` list addition.** The comment above that block
   states its entries exist because functions like `Sign` and `UnitStep` are
   "trivial comparisons" on an int64 buffer — `UnitBox` qualifies on the same
   logic. But every entry on that list is justified by a *measured* workload,
   per the file's own comment style, and `UnitBox` has none. Decision: ship
   `UnitBox` without an `AWARE` entry, and say so explicitly in the
   implementation (comment/changelog) rather than silently omitting it. This
   is a deliberate deferral, not an oversight, and can be revisited once a
   workload measurement justifies it.
3. **Closed-interval boundary confirmed.** `UnitBox[-1/2] = UnitBox[1/2] = 1`.
   This is consistent with `UnitStep[0] = 1`, which the spec doc states
   explicitly at `docs/spec/builtins/elementary-functions.md:578`. No further
   investigation needed — the user's own definition in the original request
   already specified the closed interval (`-1/2 <= x <= 1/2`).
4. **`add_test()` required for the new test target.** The suite currently has
   456 `add_executable` lines but only 221 `add_test()` lines in
   `tests/CMakeLists.txt` — roughly half the suite builds but never runs under
   `ctest`. The file's own header comment (lines 4-7) states that targets
   lacking `add_test()` "remain invisible to ctest," and `clip_tests` is one
   of the invisible ones. Decision: the new `tests/test_unitbox.c` target must
   include both `add_executable` AND `add_test()`, following the header's
   stated intent rather than repeating the `clip_tests` precedent.
