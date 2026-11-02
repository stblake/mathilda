# Plan: Implement `N[...]` in Mathilda (machine + arbitrary precision)

## Context

Mathilda currently has no numeric-evaluation function. Users can write
`1/3`, `Pi`, `Sin[1.0]`, etc., but there is no way to ask "give me a
numeric approximation of this expression." Mathematica's `N[expr]` and
`N[expr, n]` are the standard vehicle for this and are heavily used.

The numeric machinery is already *mostly* in place for machine-precision
doubles: `trig.c`, `power.c`, `logexp.c`, `hyperbolic.c`, `complex.c`
all have `get_approx`-style paths that fall through to `<math.h>` / `<complex.h>`
(`csin`, `cpow`, `clog`, `cexp`, …) when any argument is `EXPR_REAL`.
`plus.c` and `times.c` collapse to doubles when a real appears. What
is missing is (a) a function that triggers this path deliberately and
(b) numeric values for the named constants (`Pi`, `E`, `EulerGamma`, …).

Arbitrary precision (`N[expr, 50]`) has no existing support: GMP is
linked only as `mpz_t` (integers), there is no MPFR dependency, no
`EXPR_MPFR` type, and no MPFR-aware paths in any builtin. This is a
real project, not plumbing.

Outcome after this plan:

- `N[expr]` and `N[expr, MachinePrecision]` produce machine-precision
  doubles using the existing infrastructure.
- `N[expr, prec]` (for integer `prec`) uses MPFR to produce arbitrary-precision
  reals and complexes. MPFR is gated behind a `USE_MPFR=1` makefile flag
  (default on) analogous to the existing `USE_ECM=1`. When disabled,
  the 2-arg form falls back to machine precision with a warning.
- Named constants (`Pi`, `E`, `EulerGamma`, `Catalan`, `GoldenRatio`,
  `Degree`) become numericalizable.
- **Mathematica-style precision literals** are parsed:
  `` 3.98`50 `` → 50-digit precision real; `` 3.98``50 `` → 50-digit
  accuracy real; `` 3`30 `` → integer 3 re-interpreted at 30-digit
  precision. Requires MPFR; without it, the parser warns and emits a
  machine-precision `EXPR_REAL`.
- **Precision / Accuracy introspection** via `Precision[x]`,
  `Accuracy[x]`, `SetPrecision[x, n]`, `SetAccuracy[x, n]`.

---

## Phase 1 — Machine-precision `N[]` (no new build deps)

### New files

- `src/numeric.h`
- `src/numeric.c`

### Public surface

```c
/* src/numeric.h */
Expr* builtin_n(Expr* res);
void  numeric_init(void);           /* registers N + attributes + docstring */
```

### Algorithm

`builtin_n(res)` — called by evaluator after args are evaluated normally:

1. Accept 1 or 2 args. Valid prec values: absent, `MachinePrecision`
   symbol, integer. Non-integer or symbolic prec → return `NULL`
   (stay unevaluated).
2. Decide `mode`: `MACHINE` (phase 1 always) or `MPFR` (phase 2, with
   `bits = ceil(prec_digits * 3.3219281)`).
3. Call `numericalize(arg, mode, bits)`; the result is already a fully
   evaluated numeric form. `expr_free(res)` and return it.

`numericalize(e, mode, bits)` — recursive walk, matching the user's
chosen "descend + re-evaluate" semantics:

- `EXPR_INTEGER` → `expr_new_real((double)e->data.integer)`
  (phase 2: MPFR variant)
- `EXPR_BIGINT`  → `expr_new_real(mpz_get_d(e->data.bigint))`
  (phase 2: `mpfr_set_z`)
- `EXPR_REAL`    → `expr_copy(e)` (machine); (phase 2: `mpfr_set_d` with
  bits if mode=MPFR — this is a *promotion*, not a refinement; callers
  should be aware the double's bits below 53 are unknown)
- `EXPR_STRING`  → `expr_copy(e)`
- `EXPR_SYMBOL`  → table lookup (see "Constants" below). Unknown symbols
  pass through unchanged — matches Mathematica's `N[x]` → `x`.
- `EXPR_FUNCTION` with head `Rational[n,d]` → direct `(double)n / d`
- `EXPR_FUNCTION` with head `Complex[re,im]` → `make_complex(numericalize(re), numericalize(im))`
- `EXPR_FUNCTION` with head `HoldForm` / `Hold` / `Unevaluated` → return
  unchanged (preserve hold semantics).
- `EXPR_FUNCTION` otherwise → numericalize each arg, rebuild with the
  same head, then **`return evaluate(rebuilt)`** so `Plus[2.0, 3.0]`
  collapses to `5.0`, `Sin[1.0]` collapses to `0.841…`, etc. This is
  the key trick that reuses all existing numeric paths.

### Constants table

In `numeric.c`, a small static table keyed on symbol name:

| Symbol | Machine value |
|---|---|
| `Pi`          | `M_PI`                           |
| `E`           | `M_E` (or `exp(1.0)`)            |
| `EulerGamma`  | `0.5772156649015328606065120900824024310421` |
| `Catalan`     | `0.9159655941772190150546035149323841107741` |
| `GoldenRatio` | `(1.0 + sqrt(5.0)) / 2.0`        |
| `Degree`      | `M_PI / 180.0`                   |
| `Infinity`, `ComplexInfinity`, `Indeterminate`, `True`, `False` | pass through |

`I` stays structural (`Complex[0,1]` is already handled via the FUNCTION
case). Phase 2 extends this table with `mpfr_const_pi`, `mpfr_const_euler`,
`mpfr_const_catalan`; `E` is `mpfr_exp(1)`, `GoldenRatio` is computed
at requested precision.

### Attributes

`ATTR_PROTECTED | ATTR_LISTABLE`. **Not** `HoldAll` (Mathematica evaluates
`N`'s arguments normally — `N[Sin[Pi]]` first reduces to `N[0]`, then
`0.0`). `ATTR_LISTABLE` gives us `N[{1/3, 1/7}]` → `{N[1/3], N[1/7]}` for
free via `apply_listable` in `eval.c:109-156`.

### Registration

In `src/core.c`, add `#include "numeric.h"` and call `numeric_init()`
at the end of `core_init()`. `numeric_init()`:

```c
symtab_add_builtin("N", builtin_n);
symtab_get_def("N")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
symtab_set_docstring("N",
    "N[expr]\n\tGives a machine-precision numerical approximation of expr.\n"
    "N[expr, n]\n\tGives a numerical approximation to n decimal digits (requires MPFR).");
```

### Tests

New `tests/test_numeric.c`. Add to `tests/CMakeLists.txt`.

Phase 1 cases (string-compare via `assert_eval_eq`):

- `N[Pi]` → `3.14159` (tolerance: format is `%g`, so this is the exact
  printed form)
- `N[E]` → `2.71828`
- `N[1/3]` → `0.333333`
- `N[Sin[1]]` → `0.841471`
- `N[Sin[Pi]]` → `0.` (because `Sin[Pi]` evaluates to `0` before `N` sees it)
- `N[Plus[Pi, 1]]` → `4.14159`
- `N[{1/3, Pi, E}]` → `{0.333333, 3.14159, 2.71828}`
- `N[foo]` → `foo` (unknown symbol untouched)
- `N[foo[Pi, x]]` → `foo[3.14159, x]` (descend + re-evaluate semantics)
- `N[HoldForm[Pi]]` → `HoldForm[Pi]` (Hold preserved)
- `N[Complex[1, 2]]` → `1. + 2. I` (or whatever the printer emits — test
  via `FullForm` for determinism: `Complex[1., 2.]`)

Also add `?N` to `tests/test_information.c` if that pattern exists there.

### Critical files to modify (Phase 1)

- `src/numeric.c` — new
- `src/numeric.h` — new
- `src/core.c` — add include + init call (currently ends `core_init()` at
  the cluster of `*_init()` calls around the file's tail)
- `tests/test_numeric.c` — new
- `tests/CMakeLists.txt` — register new test
- `Mathilda_spec.md` — document `N[]`

No makefile changes in Phase 1; `SRC = $(wildcard src/*.c)` picks up
`numeric.c` automatically.

---

## Phase 2 — Arbitrary precision via MPFR (gated behind `USE_MPFR`)

### Build system

Mirror the `USE_ECM` pattern in `makefile`:

```make
USE_MPFR ?= 1
ifeq ($(USE_MPFR), 1)
  CFLAGS  += -DUSE_MPFR
  LDFLAGS += -lmpfr
endif
```

When `USE_MPFR=0`:
- `builtin_n` with a second argument prints a one-shot warning ("N[expr, n]:
  arbitrary precision unavailable, using machine precision") and falls
  back to the machine path.
- `EXPR_MPFR` is not defined, and all MPFR code is `#ifdef USE_MPFR`-guarded.

### New expression type

In `src/expr.h`, inside the existing tagged union — gated so phase-1
builds without MPFR still compile:

```c
typedef enum {
    EXPR_INTEGER,
    EXPR_REAL,
    EXPR_SYMBOL,
    EXPR_STRING,
    EXPR_FUNCTION,
    EXPR_BIGINT,
#ifdef USE_MPFR
    EXPR_MPFR,
#endif
} ExprType;
```

Payload (add to the union):

```c
#ifdef USE_MPFR
    mpfr_t mpfr;   /* carries its own precision in bits */
#endif
```

MPFR's `mpfr_t` tracks precision per-value, so we don't need a separate
precision field. Rounding mode is global-ish in Mathilda; default to
`MPFR_RNDN` (round to nearest, ties to even).

### Core helpers (`src/expr.c`)

```c
#ifdef USE_MPFR
Expr* expr_new_mpfr_bits(mpfr_prec_t bits);           /* zero-inited */
Expr* expr_new_mpfr_from_d(double d, mpfr_prec_t bits);
Expr* expr_new_mpfr_from_si(long v, mpfr_prec_t bits);
Expr* expr_new_mpfr_from_mpz(const mpz_t z, mpfr_prec_t bits);
#endif
```

Update `expr_free`, `expr_copy`, `expr_eq`, `expr_compare`, `expr_hash`
to handle `EXPR_MPFR`. Equality is bitwise on the mpfr payload at the
*larger* precision (round up, compare). Ordering between `EXPR_MPFR`
and other numerics goes through mpfr's comparison primitives.

### Printer (`src/print.c`)

```c
#ifdef USE_MPFR
case EXPR_MPFR: {
    /* Number of significant decimal digits for the value's precision. */
    long digits = (long)ceil(mpfr_get_prec(e->data.mpfr) / 3.3219280948873626);
    mpfr_printf("%.*Rg", (int)digits, e->data.mpfr);
    break;
}
#endif
```

Matches the `%g`-style existing double output. `InputForm` should use
`%.*Re` for round-trippable form.

### Numericalization (Phase 2 in `numeric.c`)

Extend `numericalize` to take `bits`. When `mode == MPFR`:

- Constants: `mpfr_const_pi`, `mpfr_const_euler`, `mpfr_const_catalan`.
  `E` = `mpfr_exp(one)`. `GoldenRatio` = `(1 + sqrt(5)) / 2` computed at
  `bits + 20` guard bits, then rounded to `bits`.
- `EXPR_INTEGER`/`EXPR_BIGINT` → MPFR at requested precision.
- `EXPR_REAL` at MPFR request → *promote* to MPFR with bits; value is
  exact to 53 bits and arbitrary beyond.
- Function recursion is identical, but rebuilt expression is re-evaluated
  and the MPFR-aware builtins (below) collapse it.

### MPFR-aware numeric builtins

Add MPFR branches to existing numeric-collapse paths:

| File | Function | New branch |
|---|---|---|
| `plus.c`  | numeric fold   | if any arg is MPFR → fold all via `mpfr_add` at max precision |
| `times.c` | numeric fold   | if any arg is MPFR → fold via `mpfr_mul` at max precision |
| `power.c` | `cpow` path (lines 318–375) | if base or exp is MPFR → `mpfr_pow`; complex case via re/im MPFR pair |
| `trig.c`  | each `csin`/`ccos`/… call | if `get_approx_mpfr` succeeds (new helper) → `mpfr_sin`, `mpfr_cos`, `mpfr_tan`, `mpfr_atan2`, etc. |
| `hyperbolic.c` | each `csinh`/… call | `mpfr_sinh`, `mpfr_cosh`, `mpfr_tanh`, `mpfr_asinh`, … |
| `logexp.c` | `clog`/`cexp` calls | `mpfr_log`, `mpfr_exp`; branch-cut guards preserved |
| `complex.c` | `Abs`, `Arg`, `Re`, `Im`, `Conjugate` | honor MPFR components |
| `arithmetic.c` | `expr_is_integer_like` stays the same; add `expr_is_numeric_like` that includes MPFR |

New helper (model after `get_approx`, existing file `trig.c:158-186`):

```c
#ifdef USE_MPFR
/* Like get_approx, but fills MPFR re/im and returns the precision used.
 * Precision = max(bits of any MPFR input, bits_default) for uniform result. */
bool get_approx_mpfr(const Expr* e,
                     mpfr_t re, mpfr_t im,
                     mpfr_prec_t bits,
                     bool* is_inexact);
#endif
```

Lives in `numeric.c` and is extern'd in `numeric.h` so all four trig/hyp/log
modules can use it.

### Significance arithmetic (v1)

Phase 2 v1 uses `max(prec)` of inputs for each operation's output
precision — simple, predictable, and MPFR-native. True Mathematica-style
significance arithmetic (where precision decreases through lossy ops)
is deferred; documented in `Mathilda_spec.md` as a known deviation.

### Precision literal syntax in the parser

Extend `parse_number` in `src/parse.c:88` to recognize a trailing
precision/accuracy marker after a completed numeric mantissa:

```
<number> ::= <mantissa> [ "`" [ <prec-number> ]           /* precision */
                       | "``" <acc-number>                /* accuracy  */
                       ]
```

Semantics (matching Mathematica):

| Literal   | Meaning |
|-----------|---------|
| `3.98`      | machine-precision double |
| `` 3.98`50 ``  | MPFR real, 50-digit precision (bits = ⌈50·log₂10⌉ = 167) |
| `` 3.98` ``    | machine-precision double with `MachinePrecision` tag (functionally same as `3.98`) |
| `` 3.98``50 `` | MPFR real, 50-digit *accuracy* |
| `` 3`30 ``     | integer 3 widened to MPFR at 30-digit precision |
| `` 3.0`15 ``   | MPFR at 15-digit precision (less than machine, down-rounded) |

**Disambiguation.** The backtick already means "context separator" in
identifiers (`src/parse.c:75,468,483`). There is no ambiguity because
the precision form is only consumed *inside `parse_number`*, after a
valid numeric mantissa — the identifier lexer never sees a digit first.

**USE_MPFR=0 fallback.** When MPFR is not compiled in, `parse_number`
emits a one-shot warning "precision literal ignored: USE_MPFR=0" on
the first `` ` ``-suffixed literal, and returns the mantissa as
`EXPR_REAL`. The precision suffix is consumed but discarded.

**Grammar notes.** The parser must *not* consume a trailing `` ` ``
that is immediately followed by an identifier character
(`isalpha` / `_` / `$`) — that would be a Mathematica context path,
handled by the existing path in `src/parse.c:468-483`. So the
decision is "backtick followed by digit/`.`/`+`/`-`" → precision marker;
"backtick followed by alpha/`_`/`$`" → context; anything else → error
or end-of-number.

### Precision / Accuracy builtins

New file `src/precision.c` (+ `.h`) implementing four builtins, all
registered in `numeric_init()` and given
`ATTR_PROTECTED | ATTR_LISTABLE`:

| Builtin | Behavior |
|---|---|
| `Precision[x]`        | `EXPR_INTEGER`/`EXPR_BIGINT`/`Rational` → symbol `Infinity`; `EXPR_REAL` → symbol `MachinePrecision`; `EXPR_MPFR` → decimal digits = `mpfr_get_prec(x) / log₂10`; `Complex[a,b]` → `min(Precision[a], Precision[b])`; non-numeric → `Infinity` or symbolic `Precision[x]` (Mathematica uses the former for atoms). |
| `Accuracy[x]`         | Integers/rationals → `Infinity`; `EXPR_REAL` → `15.9546 - log10(|x|)` clamped at zero for `x=0` → returns `Infinity`; `EXPR_MPFR` → `Precision[x] - log10(|x|)`; `Complex[a,b]` → `min(Accuracy[a], Accuracy[b])`. |
| `SetPrecision[x, n]`  | Recursively walks like `N`. `n = MachinePrecision` down-converts MPFR → `EXPR_REAL`; integer `n` promotes reals/ints to MPFR at `n` digits. Exact numbers (integer, rational) are promoted to MPFR when `n` is finite; kept exact only when `n=Infinity`. Pads with zero bits on promotion from `EXPR_REAL` → MPFR (matches Mathematica's "pad zeros" semantics for machine→bignum). |
| `SetAccuracy[x, n]`   | Like `SetPrecision`, but derives MPFR precision as `bits = ⌈(n + log10(|x|)) · log₂10⌉`; zero-magnitude values get `bits = ⌈n · log₂10⌉`. |

All four follow the same `numericalize`-style recursive descent as `N`,
so `SetPrecision[Sin[Pi/4] + 1/3, 30]` lifts everything to 30-digit MPFR
and re-evaluates. The shared walker lives in `numeric.c` and is
parameterized by a "conversion policy" callback (MPFR bits target,
or "down-to-machine" marker).

### Tests (Phase 2)

Extend `tests/test_numeric.c`:

- `N[Pi, 50]` → `3.1415926535897932384626433832795028841971693993751`
- `N[E, 30]` → `2.71828182845904523536028747135`
- `N[Sin[1], 40]` → `0.8414709848078965066525023216302989996226`
- `N[Sqrt[2], 100]` → 100-digit value
- `N[1/7, 30]` → 30 digits of `0.142857...`
- `N[Pi + E, 20]` → 20 digits of `5.85987448204883...`
- `N[Exp[I Pi], 30]` → `-1. + 0. I` (or close, with tolerance for last bit)
- Round-trip: `ToExpression[ToString[N[Pi, 40], InputForm]]` equals
  `N[Pi, 40]` (if `ToString[..., InputForm]` exists; otherwise skip)

Precision-literal parsing:

- `` 3.98`50 `` → `3.9800000000000000000000000000000000000000000000000`
- `` Precision[3.98`50] `` → `50.`
- `` Accuracy[3.98`50] `` → `49.4...` (approximately)
- `` 3`30 `` → `3.00000000000000000000000000000` (30 digits)
- `` Precision[3.14] `` → `MachinePrecision`
- `` Precision[1/3] `` → `Infinity`
- `` Precision[3] `` → `Infinity`

`SetPrecision` / `SetAccuracy`:

- `SetPrecision[Pi, 40]` → 40-digit Pi
- `SetPrecision[1/3, 20]` → 20-digit MPFR `0.33333...`
- `SetPrecision[3.14159, 50]` → `3.1415899...` padded with zero bits
  beyond bit 53 (matches Mathematica)
- `Precision[SetPrecision[Pi, 40]]` → `40.`
- `SetAccuracy[0.1, 30]` → `0.1` with accuracy 30 digits
- `SetPrecision[Pi, MachinePrecision]` → `EXPR_REAL` round-trip
- Listability: `SetPrecision[{Pi, E}, 30]` → `{N[Pi,30], N[E,30]}`

Because arbitrary-precision outputs are long strings, use a
helper `assert_eval_startswith` (new, in `tests/test_utils.h`)
for tolerant prefix matches.

Add a test for the fallback path: build with `USE_MPFR=0`, confirm
`N[Pi, 50]`, `` 3.98`50 `` and `SetPrecision[Pi, 50]` each produce a
warning and a machine-precision result.

### Critical files to modify (Phase 2)

- `makefile` — USE_MPFR flag
- `src/expr.h` — `EXPR_MPFR` enum + `mpfr_t` union payload (both `#ifdef`-guarded)
- `src/expr.c` — constructors, `expr_free`, `expr_copy`, `expr_eq`,
  `expr_compare`, `expr_hash`
- `src/parse.c` — extend `parse_number` (line 88) for `` ` `` / ``` `` ``` precision/accuracy markers
- `src/print.c` — `EXPR_MPFR` case
- `src/numeric.c` — MPFR numericalize path, `get_approx_mpfr`, shared
  "set precision/accuracy" walker
- `src/precision.c` / `src/precision.h` — new: `Precision`, `Accuracy`,
  `SetPrecision`, `SetAccuracy` builtins
- `src/plus.c` — MPFR fold branch
- `src/times.c` — MPFR fold branch
- `src/power.c` — MPFR branch in numeric path (currently around lines 318–375)
- `src/trig.c` — MPFR variants of sin/cos/tan/cot/sec/csc and their inverses
- `src/hyperbolic.c` — MPFR variants
- `src/logexp.c` — MPFR variants of Log, Exp, Log[b, z]
- `src/complex.c` — MPFR-aware Abs/Arg/Re/Im/Conjugate
- `src/arithmetic.c` — `expr_is_numeric_like` helper if needed
- `src/core.c` — register `MachinePrecision` as a protected symbol
- `tests/test_numeric.c` — phase-2 cases including precision literals,
  `Precision`/`Accuracy`, `SetPrecision`/`SetAccuracy`
- `tests/test_utils.h` — prefix-match helper
- `Mathilda_spec.md` — document MPFR behaviour, precision literal syntax,
  Precision/Accuracy semantics, significance-arithmetic deviation
- `README.md` — note the MPFR build dep (optional)

---

## Verification

Phase 1:

1. `make clean && make -j$(nproc)` — clean build with no warnings.
2. `./Mathilda` REPL smoke test:
   ```
   N[Pi]
   N[Sin[1]]
   N[1/3]
   N[{Pi, E, 1/7}]
   N[foo[Pi, x]]
   ```
3. `cd tests/build && cmake .. && make && ./numeric_tests` — all pass.
4. `valgrind --leak-check=full ./Mathilda <<< 'N[Sin[Pi/4] + Cos[Pi/3]]; Quit[]'`
   — zero leaks (the "definitely lost" and "indirectly lost" sections
   must both be zero).

Phase 2:

5. `USE_MPFR=1 make -j$(nproc)` — clean build.
6. REPL:
   ```
   N[Pi, 50]
   N[Exp[1], 100]
   N[Sqrt[2], 30]
   N[Sin[Pi/3], 40]
   3.98`50
   Precision[3.98`50]           (* expect 50. *)
   Accuracy[3.98`50]            (* expect ~49.4 *)
   SetPrecision[Pi, 40]
   Precision[SetPrecision[1/3, 25]]   (* expect 25. *)
   SetAccuracy[0.1, 30]
   ```
7. `USE_MPFR=0 make -j$(nproc)` — clean build; `N[Pi, 50]`,
   `` 3.98`50 ``, and `SetPrecision[Pi, 50]` each emit the fallback
   warning and return a double.
8. Run all existing tests plus `test_numeric` under both USE_MPFR
   settings — no regressions.
9. Valgrind under MPFR with several `N[…, prec]` and `SetPrecision[…]`
   calls — zero leaks.

## Out of scope (future follow-ups)

- True Mathematica-style significance arithmetic (precision degradation
  through lossy operations).
- `$MachinePrecision`, `$MaxPrecision`, `$MinPrecision` system variables.
- `Glaisher`, `Khinchin` and other less-common constants.
- `NSolve`, `NIntegrate`, `NDSolve` — all need numeric infrastructure
  but are separate projects.
