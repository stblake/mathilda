# PolynomialReduce implementation

Builtin `PolynomialReduce[poly, {p1..pn}, {x1..xk}]` → `{{a1..an}, b}` with
`Σ ai·pi + b == poly` and `b` fully reduced. Domains: Q, RationalFunctions
(Q(params)), Modulus (GF(p)). Integers / InexactNumbers deferred (warn+decline).

## Tasks
- [x] Shared named-order → weight-matrix helper in `groebner.{c,h}`
      (`gb_build_order_matrix`, `gb_classify_named_order`); refactored
      `monomials.c` to use it.
- [x] Fixed `GroebnerBasis` MonomialOrder parser (DegreeLexicographic now
      honoured via the shared helper; Negative* → validated → Lex fallback).
- [x] `gfp_divmod` (cofactor-tracking) in `gbmod.{c,h}`.
- [x] `flint_polynomial_reduce` in `flint_bridge.{c,h}` (fmpq_mpoly_divrem_ideal).
- [x] `src/poly/polynomialreduce.c`: front-end, Q path (FLINT Lex + gb_divmod),
      RationalFunctions Q(params) engine, Modulus path, deferred-domain decline.
- [x] Wire: `poly.h` decls, `poly_init()` call, `SYM_PolynomialReduce`,
      docstring in `info.c`, default `Options` in `options_builtin.c`,
      `ATTR_PROTECTED`.
- [x] `tests/test_polynomialreduce.c` + CMake registration + COMMON_SRC entry.
- [x] Docs: `structural-manipulation.md` section + changelog entry.
- [x] Build clean (`make`), tests pass, spec examples replayed, valgrind-clean
      (matches baseline), graph self-review.

## Review

**What shipped.** `PolynomialReduce[poly, {pi}, {xi}] -> {{ai}, b}` with
`Σ ai·pi + b == poly` and `b` reduced. Domains: Rationals over Q,
RationalFunctions over Q(params) (default), Modulus (GF(p)). Integers /
InexactNumbers / nonzero Tolerance are deferred (accepted, then decline with a
note). Plus a `GroebnerBasis` fix: `DegreeLexicographic` (and the shared named
orders) are now honoured instead of silently downgrading to Lex.

**Engines.** gb_divmod (`groebner.c`) is the authoritative Q engine; the pure-Q
Lexicographic case uses FLINT `fmpq_mpoly_divrem_ideal` (verified to match
gb_divmod and Wolfram). RationalFunctions uses an Expr-coefficient division over
Q(params) with `Together`/`Cancel` field arithmetic (reuses `CoefficientRules`).
Modulus uses the new `gfp_divmod`.

**Verification.** Every applicable spec example reproduced (pure-Q and DegRevLex
match Wolfram *exactly*; RationalFunctions remainders match, quotients match up
to `GroebnerBasis`'s pre-existing sign normalisation; the reconstruction
invariant `q.p + b == poly` holds in every domain). Fixed one real leak
(`field_neg` was borrowing an owned `field_mul` temporary); valgrind now matches
the leak baseline exactly. `polynomialreduce_tests`, `groebner_tests`,
`coefficient_rules_tests` all pass.

**Deferred / known limitations.** Integers & InexactNumbers+Tolerance domains
(distinct ring/numeric engines). FLINT fast path enabled for Lexicographic only
(DegLex/DegRevLex use gb_divmod pending grevlex-convention verification on a
FLINT build). Modulus not combined with parameters.

---

# BitLength implementation

Builtin `BitLength[n]` -> number of binary bits to represent integer `n`.
`n>0`: `Floor[Log[2,n]]+1`; `n=0`: `0`; `n<0`: `BitLength[BitNot[n]]`
(`BitNot[n]=-n-1`). Attributes `Listable`, `Protected`; arity 1. Seeds a new
`src/bitwise/` module.

## Tasks
- [x] `src/bitwise/{bitwise.h,bitwise.c,bitlength.c}` — hub + GMP builtin +
      `::int`/`::argx` diagnostics.
- [x] Build wiring: `makefile` (SRC wildcard + `-I./src/bitwise`), `core.c`
      (`#include` + `bitwise_init()`).
- [x] Interned symbol `SYM_BitLength` (`sym_names.{h,c}`); attributes in
      `attr.c`; docstring in `info.c`.
- [x] Packed/NDArray: `ndk_BitLength_ii` int64 kernel + descriptor +
      registration (`ndinteger.c`); `AWARE` + `INT64_OK` (`pack.c`).
- [x] Compile[]: opcode `OP_BLEN_I` (`compile_internal.h`), `int_only_head`
      (`compile_infer.c`), codegen arm (`compile_emit_int.c`), VM helper
      `int_blen` + dispatch (`compile_vm.c`). Updated the "Bit* absent" notes
      (`compile_internal.h`, `docs/design/compile_state.md`).
- [x] Tests: `tests/test_bitwise.c` (29 cases) + CMake target; leak script
      `tests/scripts/bitwise_leakcheck.sh`.
- [x] Docs: new `docs/spec/builtins/bitwise.md`, `Mathilda_spec.md` overview
      row, `docs/spec/changelog/2026-09-14.md` note.

## Review
**Verification.** REPL: positives/negatives/bignums/INT64_MIN all correct;
`BitLength[n]==IntegerLength[n,2]` for n=1..500; packed `Range` and visible
`NDArray[..int64..]` agree and are single-headed; `Compile` scalar & rank-1
array both `Compiled -> True`; `100!` -> 525. `bitwise_tests` (29) pass;
`bitwise_leakcheck.sh` reports 0 leaks under `leaks`. Audits green:
`check-c99`, `check-packed-aware`, `check-array-exactness`, `check-nd-surfaces`.
`check-compile-coverage` fails only on the pre-existing image/fit baseline
(documented RED-on-main; `BitLength` is NOT in the flagged list).

**Notes.** ND/Compile/pack code stays in its canonical homes (`ndinteger.c`,
`src/compile/`, `pack.c`) though the scalar builtin lives in `src/bitwise/` —
mirrors how `numbertheory/` heads are wired. int64 kernel complements the
unsigned magnitude, so `INT64_MIN` is exact with no overflow (cleaner than
`IntegerLength`'s negate path).
