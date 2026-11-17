# Numeric-tower review (2026-11): machine, bignum and arbitrary precision

A review of the three number layers and the boundaries between them:
machine integers (`EXPR_INTEGER`, `int64_t`) and GMP bignums (`EXPR_BIGINT`),
exact rationals, machine reals (`EXPR_REAL`, IEEE double), and MPFR
arbitrary-precision reals (`EXPR_MPFR`). It looked for crashes, wrong
answers and inconsistencies between surfaces (scalar, packed, `NDArray`,
`Compile[]`).

Method: three code-mapping passes over `src/` and about 190 probes against
the binary. Expected values come from first principles: exact identities,
the IEEE rules, and GMP/MPFR at higher precision. The regression corpus is
`tests/test_numeric_tower_review.c` (ctest `numeric_tower_review_tests`).
Each test there is keyed by the finding ids below.

**Scope note.** Mathilda's arbitrary-precision numbers are deliberately
**fixed precision**: the MPFR bit width is the precision, and no error bound
is carried. The project does not use significance arithmetic. A
fixed-precision result such as `(1 + 10^-30) - 1.0`40`, which keeps its
40-digit label after the cancellation, is therefore expected behaviour, not a
defect. That candidate finding (D1) was withdrawn.

## The machine/software boundary policy (v0.6594)

> **A machine number has a 53-bit mantissa and a normal exponent.**
>
> A machine-arithmetic result outside `[DBL_MIN, DBL_MAX]` is the same value
> at machine precision: a 53-bit MPFR with MPFR's exponent range. That covers
> overflow, an underflow to zero from nonzero operands, and a subnormal.
>
> `Overflow[]` and `Underflow[]` are reserved for values past the MPFR range.

`Power` and `Exp` already followed this rule: `Exp[1000.]` is
`1.97*10^434`. Everything else did not. `Plus` and `Times` returned
`Overflow[]` or a silent `0.`. The parser produced a raw `inf.0`. `N[]` and
the parser produced subnormal doubles with a few good digits.

The rule is now applied at every choke point:

| Choke point | Where |
|---|---|
| machine `Plus` / `Times` | `add_numbers` (plus.c), `multiply_numbers` (times.c) → `numeric_machine_binop_promote` |
| machine call results | `numeric_promote_result_if_degenerate`, `numeric_machine_real_or_promote` |
| promoted value demoted back | `demote_to_machine_if_fits` (it demoted subnormals, contrary to its comment) |
| N of an exact value | `mpfr_to_machine`, the rational leaf in `numericalize` |
| machine literals | `machine_real_from_decimal` (parse.c) |

The predicate is `numeric_machine_subnormal` (numeric.h).

## Findings

Status: **fixed** (version), **withdrawn** (with the reason), or **open**.

### A. Crashes and undefined behaviour
| # | Repro | Was | Status |
|---|---|---|---|
| A1 | `Sqrt[-2^63]` | segfault (`-INT64_MIN` in power.c) | fixed: `(-1)^r 2^(63 r)` |
| A2 | `Quotient[1.*^300, 1]` | `-2^63` (UB cast) | fixed: exact via `mpz_set_d` |
| A3 | `N[1, 10^9]` | 8 s timeout | withdrawn: a 10^9-digit result is that slow to print; not a hang |
| A4 | `Compile` `Norm[v, 1.*^300]` | UB `(int64_t)` cast | fixed: stays Real |
| A5 | `3^(2^62)`, `(2^1000/3)^(10^9)` | `mpz_pow_ui` aborts the process | fixed: `General::ovfl`, `Overflow[]` past GMP's range |
| A6 | `x^(-2^63/3)` | the printer recursed until the stack overflowed | fixed: `make_rational_neg` |
| — | `Times` fusion pass | `INT64_MIN % -1` traps (SIGFPE) | fixed: skipped |

### B. Exact layer
| # | Repro | Was | Status |
|---|---|---|---|
| B1 | `Divide[1/3^39, 3^39]` | `-1/7710839681495117703` (int64 wrap) | fixed |
| B2 | `GCD[-2^63]`, `LCM[-2^63, 3]` | negative | fixed: GMP fold; `lcm()` saturates instead of wrapping |
| B3 | `(-3)^-41` | `Rational[1, -3^41]` | fixed: `bigint_pow` builds the canonical reciprocal |
| B4 | `Rational[2^70, 2^70]` | not reduced | fixed |
| B5 | `000000000000000000001` | a BigInt; `Hash`, `Union` and `Association` disagreed with `SameQ` | fixed: parser normalises; an int64-sized BigInt hashes as Integer |
| B6 | `(-1)^(2^70)` | unevaluated | fixed: parity |
| B7 | `get_numeric_value` (sort.c) | `mpz_get_d(n)/mpz_get_d(d)` = NaN | fixed: `numeric_machine_value` |

`expr_bigint_normalize` now uses a portable int64 test
(`expr_mpz_get_int64`). `mpz_fits_slong_p` tests a 32-bit `long` on LLP64
Windows. The other ~87 `mpz_fits_slong_p` sites are only conservative there,
not wrong. They are left as they are.

### C. Machine ↔ software
| # | Repro | Was | Status |
|---|---|---|---|
| C1 | `1.*^400`, `1.*^-400`, `1.*^-320` | `inf.0`, `0.`, a 4-digit subnormal | fixed (policy) |
| C2 | `1.*^308 + 1.*^308`, `1.*^-200 * 1.*^-200` | `Overflow[]`, `0.` | fixed (policy) |
| C3 | packed `{1.*^308} * 10` | `{inf.0}` | open: see section E |
| C4 | `Sqrt[-1.]` (all surfaces) | `6.1*^-17 + 1. I` | fixed: `mth_sincospi`, shared by scalar and ND |
| C5 | `Log[0.]` | `Indeterminate` | withdrawn: no ground truth available. `Indeterminate` is defensible for an approximate zero |
| C6 | `N[Pi,30]^1.5`, `Complex[1.5`30, 2.]` | 30 digits kept | fixed: contagion in `Power`; `Complex[]` is one kind |
| C7 | `1.5 < 10^400`; `Sort[{N[10^500], N[10^400]}]` | unevaluated; wrong order | fixed: exact mixed-kind order |
| C8 | `N[10^-320]`, `2.^-1074` | subnormal machine doubles | fixed (policy) |
| C9 | `1.*^15` prints `1e+15` | OutputForm has no decimal point | open: see section E |
| C10 | int64 overflow in packed / `Compile` | — | verified sound (`checked_int.h`) |

### D. Arbitrary precision (fixed-precision model)
| # | Repro | Got | Status |
|---|---|---|---|
| D1 | no significance arithmetic | — | withdrawn: design choice |
| D2 | `Precision[0``20]` | `20.` | withdrawn: under the fixed model a zero's precision is its MPFR width |
| D3 | `Accuracy[N[10^-400, 20]]`; `Accuracy[N[10^500, 20]]` | `20.` (should be 420); `Overflow[]` | open: `expr_log10_abs` uses `mpfr_get_d`/`mpz_get_d` |
| D4 | `Precision[SetAccuracy[Pi 10^50, 10]]` | `10.` (should be ~60) | open: a symbolic value has magnitude 0 |
| D5 | `Exp[-N[10^10, 30]]` (below MPFR's range) | `0.e-30`, no message | open (decided 2026-11-17: `Underflow[]` + `General::unfl`) |
| D6 | `N[Exp[10^9], 20]` | `Overflow[]` | decided 2026-11-17: keep the library default; no change |
| D7 | `N[1/3, -1]`, `SetPrecision[1.5, 0]`, `$MinPrecision` | unevaluated, no message | open |
| D8 | `` 1.5`+nan ``; ``` 0``-5 ``` | `2.` (accepted silently); accuracy clamped to 1 | open: parser `strtod` in the precision mark |
| D9 | `-2.`30^2`, `-2`10^2`, `-2.5*^3^2` | `4.`, `4.`, `6.25*10^6` (should be negative) | open: **wrong answer**. `number_followed_by_power` does not skip a precision mark or `*^` |
| D9 note | `Exp[-10.`30^10]` | `Overflow[]` | the same parser bug: it read as `Exp[(-10.`30)^10]`. It was recorded as an underflow finding at first |

### E. Found while fixing (2026-11-17)
| # | Repro | Got | Status |
|---|---|---|---|
| C3 | packed / NDArray `{1.*^308, 2.} * 10`, `+`, `{1.*^-200} * 1.*^-200`, `Total` | `{inf.0, 20.}`, `{0., ...}`, `Overflow[]` | open: **wrong answer**. `nd_ew_chunk` has no degenerate-lane decline; the scalar path is fixed |
| E5 | packed `{1., 2.} * 10^400`, `+ N[10^400, 30]` | `{inf.0, inf.0}` | open: `real_scalar_value` (ndarray.c) decodes with `mpz_get_d`/`mpfr_get_d` |
| E9 | `Compile[{{x, _Real}}, x x][1.*^-200]` | `0.` (the interpreter gives `1.*10^-400`) | open: Compile scalar/rank-1 underflow; overflow already falls back |
| E4 | `Quotient[N[10^30, 40] + 1/2, 1] - 10^30` | `19884624838656` (should be 0) | open: **wrong answer**. The `Quotient` real fallback reads MPFR through a double |
| E2 | `Exp[1000.] * 0.` | `0.e-15` (an MPFR zero) | open: an MPFR value times a machine `0.` should be a machine `0.` |
| E1 | `(1.*^308 * 10.)/1.*^10` | 53-bit MPFR, `MachineNumberQ` False | open (decided 2026-11-17: demote once in range) |
| E3 | `(2/3)^(-2^63)` | `1/(2/3)^9223372036854775808` | open, cosmetic: should be `(3/2)^(2^63)` |
| C9 | `1.*^15`, `$MachineEpsilon` | `1e+15`, `2.22045e-16`; `Precision[2.^-1074]` is `15.9546` | open (decided 2026-11-17: always print a decimal point; a promoted value reports `MachinePrecision`) |

### F. Hygiene (no repro: from the code read)
| # | What | Status |
|---|---|---|
| F1 | about 149 `mpz_get_d` calls truncate instead of rounding to nearest (`numeric_machine_value` rounds correctly) | open: audit the sites that reach values users see |
| F2 | about 87 `mpz_fits_slong_p` / `mpz_get_si` sites test a 32-bit `long` on LLP64 Windows | open: audit; `expr_mpz_get_int64` is the portable form |
| F3 | `cpow` cleanup drops an imaginary part below 1e-14 (`(1. + 1.*^-15 I)^0.5`) | kept: within the machine error |
| F4 | `checked_int.h` cites `tests/test_compile_overflow.c`, which does not exist | open: write it or drop the reference |
| F5 | no automated check guards the machine boundary across surfaces | open: an `ND`/Compile boundary sweep (see Phase 8 in `tasks/todo.md`) |

## Known, not changed
- A promoted value does not demote back when a later result fits again:
  `(1.*^308 * 10.)/1.*^10` stays a 53-bit MPFR `1.*10^299`.
- `(1. + 1.*^-15 I)^0.5` drops its 5*10^-16 imaginary part (cpow cleanup).
  That part is within the machine error of `cpow`.
