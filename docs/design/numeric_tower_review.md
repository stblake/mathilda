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
| C3 | packed `{1.*^308} * 10` | `{inf.0}` | fixed v0.6595 (see section E) |
| C4 | `Sqrt[-1.]` (all surfaces) | `6.1*^-17 + 1. I` | fixed: `mth_sincospi`, shared by scalar and ND |
| C5 | `Log[0.]` | `Indeterminate` | withdrawn: no ground truth available. `Indeterminate` is defensible for an approximate zero |
| C6 | `N[Pi,30]^1.5`, `Complex[1.5`30, 2.]` | 30 digits kept | fixed: contagion in `Power`; `Complex[]` is one kind |
| C7 | `1.5 < 10^400`; `Sort[{N[10^500], N[10^400]}]` | unevaluated; wrong order | fixed: exact mixed-kind order |
| C8 | `N[10^-320]`, `2.^-1074` | subnormal machine doubles | fixed (policy) |
| C9 | `1.*^15` prints `1e+15` | OutputForm has no decimal point | **fixed v0.6604**: see section E (C9) |
| C10 | int64 overflow in packed / `Compile` | — | verified sound (`checked_int.h`) |

### D. Arbitrary precision (fixed-precision model)
| # | Repro | Got | Status |
|---|---|---|---|
| D1 | no significance arithmetic | — | withdrawn: design choice |
| D2 | `Precision[0``20]` | `20.` | withdrawn: under the fixed model a zero's precision is its MPFR width |
| D3 | `Accuracy[N[10^-400, 20]]`; `Accuracy[N[10^500, 20]]` | `20.` (should be 420); `Overflow[]` | fixed v0.6596: `*_get_d_2exp` |
| D4 | `Precision[SetAccuracy[Pi 10^50, 10]]` | `10.` (should be ~60) | fixed v0.6596: machine-value magnitude; digits rounded up; `0``a` when no digit is left |
| D5 | `Exp[-N[10^10, 30]]` (below MPFR's range) | `0.e-30`, no message | **fixed v0.6604**: `Underflow[]` (`-Underflow[]` for a negative zero) with `General::unfl`. The evaluator scopes MPFR's underflow flag to each `NumericFunction` call (saved, cleared, restored) and converts a zero result raised under it (`numeric_underflow_if_zero`); the machine path's MPFR recompute decides for itself (`numeric_unfl_deferred`), so the message fires once. A true zero stays a zero |
| D6 | `N[Exp[10^9], 20]` | `Overflow[]` | decided 2026-11-17: keep the library default; no change |
| D7 | `N[1/3, -1]`, `SetPrecision[1.5, 0]`, `$MinPrecision` | unevaluated, no message | fixed v0.6596: `N::precbd` / `SetPrecision::precbd`; `$MaxExtraPrecision` Infinity / Rational honoured. `$MinPrecision`/`$MaxPrecision`: D7b, implemented v0.6597 (constructor clamp, bounded requests, validated assignment) |
| D8 | `` 1.5`+nan ``; ``` 0``-5 ``` | `2.` (accepted silently); accuracy clamped to 1 | `nan`: fixed v0.6596 (strict decimal scan). ``` 0``-5 ```: withdrawn -- a zero's accuracy is its MPFR width (>= 1 bit), so a negative accuracy has no representation in the fixed model |
| D9 | `-2.`30^2`, `-2`10^2`, `-2.5*^3^2` | `4.`, `4.`, `6.25*10^6` (should be negative) | fixed v0.6595 |
| D9 note | `Exp[-10.`30^10]` | `Overflow[]` | the same parser bug: it read as `Exp[(-10.`30)^10]`. It was recorded as an underflow finding at first |

### E. Found while fixing (2026-11-17)
| # | Repro | Got | Status |
|---|---|---|---|
| C3 | packed / NDArray `{1.*^308, 2.} * 10`, `+`, `{1.*^-200} * 1.*^-200`, `Total` | `{inf.0, 20.}`, `{0., ...}`, `Overflow[]` | fixed v0.6595: degenerate-lane decline + delist (elementwise), `nd_guard` (reductions), stats fast paths, `Dot` |
| E5 | packed `{1., 2.} * 10^400`, `+ N[10^400, 30]` | `{inf.0, inf.0}` | fixed v0.6595: correctly rounded decode; an unrepresentable scalar declines |
| E9 | `Compile[{{x, _Real}}, x x][1.*^-200]` | `0.` (the interpreter gives `1.*10^-400`) | fixed v0.6595: underflow flag / subnormal result falls back (CompiledFunction, Table) |
| E12 | `Compile[{{x, _Real}, {y, _Real}}, x^y][0.5, 10^400]`; `Compile[{{x, _Real}}, 2. x][1/10^400]` | `0.` (interpreter `Underflow[]`), `0.` (interpreter `2.*^-400`) | **fixed v0.6610** (found by `check-machine-boundary`): a `_Real` argument with no machine representation (`inf.`, or `0.`/subnormal from a nonzero value) declines the compiled run (`cf_to_double`), so the interpreter keeps the exact value. E9's flag could not see it: C99 `pow(0.5, inf)` is an exact `+0` with no underflow exception |
| E4 | `Quotient[N[10^30, 40] + 1/2, 1] - 10^30` | `19884624838656` (should be 0) | fixed v0.6595: `Floor[(m - d)/n]` at the operands' precision |
| E2 | `Exp[1000.] * 0.` | `0.e-15` (an MPFR zero) | fixed v0.6595 |
| E1 | `(1.*^308 * 10.)/1.*^10` | 53-bit MPFR, `MachineNumberQ` False | **fixed v0.6604**: every `NumericFunction` result (Plus and Times included) that is a 53-bit MPFR, or a Complex of 53-bit parts, and fits a normal double is the double (`numeric_demote_machine_precision`, in the evaluator net) |
| E10 | `Range[1.*^-200, 2.*^-200, 1.*^-200]`; `Table` over it | `Range::toobig` (~10^14 elements); `Table` hung | fixed v0.6595: the iterator's bound slack is relative (`iter_real_slack`) |
| E11 | `Mean`/`Variance`/`RootMeanSquare`/`Dot` of `{1.*^308, 1.*^308}` (plain lists too) | `Overflow[]`; underflow gave `0.` | fixed v0.6595: each machine fast path declines a degenerate result |
| E3 | `(2/3)^(-2^63)` | `1/(2/3)^9223372036854775808`; `x^(-2^64)` printed `x^-18446744073709551616` in InputForm and was not put in a denominator | fixed (phase 7): past the exponent cap, a negative power of a Rational moves to the reciprocal base (`(3/2)^9223372036854775808`); a negative BigInt exponent prints the way a machine-integer one does (InputForm `x^(-n)`, OutputForm `1/x^n`, TeX `\frac`) |
| C9 | `1.*^15`, `$MachineEpsilon` | `1e+15`, `2.22045e-16`; `Precision[2.^-1074]` is `15.9546` | **fixed v0.6604**: OutputForm `1.0*^15`, `2.22045*^-16` (six digits, always a decimal point, the input syntax `*^`; arbitrary-precision display uses `*^` too; TeX `\times 10^{e}`). A 53-bit MPFR reports `MachinePrecision` and prints as a machine real: six digits for display, shortest round-trip digits and no precision mark in InputForm, which re-reads to the same value |

### S. Found by the stress corpus (2026-11-17, `tests/test_numeric_tower_stress.c`)
318 boolean cases whose truth follows from exact mathematics; 309 pass. The 9
failures are these open bugs, marked XFAIL in the suite by id:

| # | Repro | Got | Status |
|---|---|---|---|
| S1 | `I^(2^64 + 1)` | unevaluated | fixed v0.6598: period-4 arm |
| S2 | `Exp[I Pi 1.] == -1` | `False` (an imaginary part of 1.2*^-16 against an exact 0) | **accepted as is** (user, 2026-11-17): complex `Equal` compares parts separately, not \|a - b\| relative to \|a\| |
| S3 | `(1.*^200 + 1.*^200 I)^2`, `Exp[1000. + 0.5 I]` | `Overflow[]` | fixed v0.6598: the evaluator net + complex degeneracy in the promoter; integer powers by binary powering in MPFR |
| S4 | packed `{1.*^-200}^2` | `{0.}` | fixed v0.6598: fused lane check in the power kernels; ndarray_take_degenerate -> delist |
| S5 | `Exp[packed {-800., 0.}]` | `{0., 1.}` | fixed v0.6598: fused check for decaying heads (an opt-in list, replaced in v0.6612 by default-on lane modes: G5/G6); the evaluator net for every scalar fast path |
| S6 | `Product[1.*^-100, {i, 5}]` | `0.` | fixed v0.6598: compiled Sum/Product read the IEEE underflow flag |
| S7 | `16^^FF`, `2^^1010` | syntax error | fixed v0.6598: b^^digits[.digits] in parse_number |

### K. Found by the Kahan-style corpus (2026-11-17, `tests/test_numeric_kahan.c`)
307 cases: Paranoia-style machine properties, correctly rounded + - * / sqrt
against the exact result, elementary functions at hard arguments (within the
platform libm's 2-ulp bound), decimal<->binary conversion, the cancellation
classics through N, branch cuts and robust complex arithmetic, stable
statistics. 300 pass; 7 XFAIL are these two bugs:

| # | Repro | Got | Status |
|---|---|---|---|
| K1 | `-0.3/0.1`; `x/x` at `x = 1.*^-300` | `-3.0` (0.625 ulp; `Divide[-0.3, 0.1]` is the correctly rounded `-2.9999999999999996`); not exactly `1.` | **fixed v0.6600**: Times evaluates the base of `Power[b, -1]` first and fuses a machine divisor with a machine numerator on its left into `Divide` (eval.c `k1_*`); `Divide` gained a complex branch (`mth_cdiv_machine`) and a buffer kernel (`ndarray_divide`); Compile emits DIV / the buffer quotient for `x/y`. Was: `a/b` parses to `Times[a, Power[b, -1]]` and the reciprocal is rounded first, so a machine quotient is rounded TWICE. IEEE division (and x/x = 1) requires one rounding |
| K2 | `Sqrt[1.*^-200 (1 + I)]`, `Sqrt[1.*^200 - 3.*^200 I]` | relative error 5*10^-14 (about 220 ulps) | **fixed v0.6600**: `z^(1/2)` for a machine complex `z` is libm `csqrt` (`mth_csqrt_machine`) on the scalar path, in `cpow_pair` and in Compile's SQRT_C, so the surfaces bit-match. Was: exp(log(z)/2) (`cpow`), whose error grows with \|log\|z\|\| |
| K3 | shootout asinh-1-2 (`Method -> "ParallelMixedTower"`) stopped integrating at v0.6598 | `PossibleZeroQ[D[F] - f, Assumptions -> Element[x, Complexes]]` False for a correct F | **fixed v0.6604**: the cause is S3: a complex result outside the double range used to be `Overflow[]`, and the zero test re-drew such a sample. Now it is a 53-bit value whose cancellation noise read as decisive. When the machine rung makes an extended-range COMPLEX value (counted by `numeric_note_extended` at every promotion and in the evaluator net), a sample point is re-drawn again. Real extended-range values (since v0.6594) and exact inputs keep the ordinary path, which decides them (`x Pi/3^700`, `10^500 (Pi + 10^-30) - 10^500 Pi`). Broader rules tried first (re-draw on any extended value; climb one rung; climb the whole ladder) each broke one of these |

Test data found wrong on the way: `tests/test_elliptic_theta.c` expects
`EllipticTheta[1, 0.4, 0.999]` = -1.32*10^-504. Jacobi's imaginary transformation
gives |theta_1| ~ 10^-593.5, and Mathilda's arbitrary-precision value is
5.38*10^-594, stable from 30 to 120 digits. The old flushed `0.` passed the
absolute tolerance against both. The other |q| = 0.999 rows of that table are
suspect too and should be re-derived.

Not bugs (verified, the cases were adjusted): `Tan[1.*^300]` is 1.34 ulp off,
which is the platform libm's `tan`, not Mathilda; the corpus therefore uses
the libm bound of 2 ulps. `Exp[Log[z]]` at \|z\| ~ 10^+-200 and
`Tan[ArcTan[100.]]` are limited by their conditioning.

### G. Found by `check-machine-boundary` (2026-11-18, `tools/check_machine_boundary.py`)
The gate probes 161 element-wise heads at the machine-range edges over scalar /
List / packed / NDArray / Compile, and requires each lane to equal the scalar.
It found E12 (fixed v0.6610) and these, all in its BASELINE. Several are on
ORDINARY data: the second lane of every probe holds `0.5`.
| # | Repro | Got | Status |
|---|---|---|---|
| G1 | `Compile[{{x, _Integer, 1}, {y, _Integer}}, Subtract[x, y]][{10, 3}, 2]` | `{-9223372036854775796, ...}` (`x - y` is right: `{8, 1}`) | **fixed v0.6611**: integer elements negate with an integer -1 (`OP_NEG_I` / `* -1`), as `a - b` does; `OP_NEG_R` had flipped the int64 sign bit, and a real -1 turned the result into a real array (losing `2^53 + 1`) |
| G2 | `InverseJacobiDN[ToPackedArray[{0.3, 0.5}], 1.*^-10]` | real values; the scalar is `1.5708 + 12.35 I` (NDArray too) | **fixed v0.6611**: the kernel was registered real-closed, but its helper returns the complex branch `K + i y` for `0 <= v < Sqrt[1 - m]`; now escaping (computes complex, narrows only when every element is real) |
| G3 | `Mod[ToPackedArray[{0.3, 0.5}], 1.*^-10]` (NDArray, Compile rank 1) | `{0., 0.}`; the scalar is the exact fmod `9.9999978*^-11` | **fixed v0.6611** (Quotient: exact floor too, v0.6613): `m - n Floor[m/n]` floors a ROUNDED quotient (`0.5/1.*^-10` rounds to exactly `5*^9`); every surface now calls `mth_mod_machine` (exact `fmod`, moved to the divisor's sign), bit-exact on 4007 random pairs against the exact dyadic remainder. The scalar had looked right only because the evaluator's degenerate-result net recomputed its `0.`. `Quotient` keeps `Floor[m/n]` of the rounded quotient (so `Quotient[0.5, 1.*^-10]` is `5000000000`, not `4999999999`): open question |
| G4 | `ArcTan[ToPackedArray[{0.3, 0.5}]]`; ArcCot, Log2/Log10 (imaginary part), InverseErf; Compile rank-1 ArcSin/ArcCos | `0.46364760900080615` for `ArcTan[0.5]` = `0.46364760900080611621...` (scalar correctly rounded) | **fixed v0.6613**: one machine helper per head on every surface -- `src/mth_machine_elem.h` (Exp, the circular / hyperbolic functions and inverses: a real argument in the real domain takes the real libm function, else the C99 complex one with the scalar's cut side), `mth_logb_fixed_machine_z` (Log2/Log10 negative and complex), `sf_machine_binomial` for the scalar Binomial (made accurate: small-index product, tgamma quotient, lgamma only under a 1e-13 error bound -- the kernel had kept ~7 digits at `Binomial[1.*^8, 2.]`). Neither libm family is uniformly closer (`atan(0.5)` 1 ulp off, `casin(0.5)` 1 ulp off). Two consumers had depended on the old rounding: Table's exact kernel (a copy of the builtins' complex formula) and PossibleZeroQ's huge-scale gate (a 12-bit gate that a libm path could move a residual across; now the 64-bit noise level, and a rung-1 plateau certifies only a value). The gate now also probes ordinary values; it found Sinh/Tanh/Coth/Csch, Log2, Binomial and ExpIntegralE (whose kernel dropped the imaginary part for `x < 0`: registered escaping now) |
| G5 | packed `ArcTanh[{1.*^308, 0.5}]`; packed `JacobiAmplitude[2.2250738585072014*^-308, {0.5, ...}]` | real part `1.0*^-308` (subnormal); `2.2250738585072004*^-308` (subnormal) | **fixed v0.6612** with G6 (the subnormal lane declines; the List path gives the scalar's value) |
| G6 | Compile `Csch[-1.*^308]`, `Sech`, `Erfc[1.*^308]`; packed `FresnelS[{-1.*^-200, 0.5}]`, `Zeta`/`HurwitzZeta[{1.*^308, ..}, 1.*^308]`; `SphericalBesselJ[1.*^308, 10.]` on every array surface | `0.` (scalar `Underflow[]`; FresnelS `-5.24*^-601`) | **fixed v0.6612**: lane checks are the DEFAULT for every unary and binary kernel (`ND_LANE_FULL`: non-finite, subnormal, or zero from all-nonzero operands declines, as the scalar net recomputes), with two opt-outs, each a claim: `ND_LANE_NONE` (proven non-decaying: Sin, Log, ... -- the ~20% check cost on Sin) and `ND_LANE_RANGE` (exact zeros are values: Mod, Log[b, 1.], LogGamma[1.]). Compile: the kernel opcodes carry the mode in `Instr.flags` and raise `FE_UNDERFLOW`; the run check also tests `FE_OVERFLOW` (1/Sinh[800.]); the Erfc opcodes check their own lane. Values inside NIntegrate & co are left as computed. Speed within noise (Sin, Exp, Erfc, Mod, BesselJ, ArcTan2, compiled Erfc / Csch) |
| G7 | `Arg[-0.]`, `ArcCot[-0.]` | scalar `0`, `Pi/2`; kernels `Pi`, Compile `-Pi/2` | **fixed v0.6613** (decided: IEEE signed zero everywhere): `get_approx` built `x + 0.0*I`, and `-0. + 0. = +0.` lost the sign; `mth_cmplx` writes the parts directly (C99 6.2.5p13); the scalar Arg of an inexact zero goes to atan2 (`Arg[-0.] = Pi`, `Arg[0.] = 0.`, `Arg[0] = 0`) |
| G8 | `Binomial[2., -1.*^308]` | scalar `Indeterminate` (F6), kernels `0.` | **fixed v0.6613** (decided: 0.): the scalar's lattice step stops at 2^53 and `tgamma(-1.*^308 + 1)` hit its pole; it now calls the kernel's `sf_machine_binomial` |
| G9 | `Log[1.*^308, 1.*^308]` | exact `1` for machine arguments; every other surface `1.` | **fixed v0.6613** (decided: 1.): the `Log[b, b] -> 1` rule is for an exact or symbolic `b` only |
| G10 | `StruveH[{-1.*^308, 0.5}, 10.]`, Compile; `HermiteH[-1.*^308, 10.]` | over 10 s (the scalar gives up unevaluated in 7.6 s); HermiteH takes 13.5 s and returns exact `0` | **fixed v0.6613**: undefined behaviour -- seven out-of-range double -> int conversions in the machine kernels (`struve_series` ran ~2^31 terms: 7.6 s -> 1 ms; also Anger-Weber, MeijerG, the Bessel table, BesselK, Mathieu) and the four Bessel MPFR series caps `(unsigned long)(4 rhd)` (BesselK[1.*^308, 1.*^308] 18 s -> 0.46 s, NaN = decline); `hgpfq_mpfr_sum` decides an unreachable 10^6-term cap up front. And a wrong answer under it: HermiteH's 1F1 difference cancelled ~2 Sqrt(abs(n) x^2) log2(e) bits (`HermiteH[-100., 10.]` six digits, `[-1000., 10.]` `-1.0*^-1358` for `2.57*^-1608`, `[-1.*^10, 10.]` an exact 0); real x > 0 now uses `2^n U(-n/2, 1/2, x^2)`, checked against the integral representation |

### F. Hygiene (no repro: from the code read)
| # | What | Status |
|---|---|---|
| F1 | about 149 `mpz_get_d` calls truncate instead of rounding to nearest (`numeric_machine_value` rounds correctly) | fixed (phase 7): 82 calls that reach a user-visible value or decide an answer now round once to nearest (`numeric_mpz_to_double`, round-half-even in GMP integers so a no-MPFR build agrees; `numeric_mpz_ratio_to_double` for n/d pairs, which read inf/inf = NaN for a huge Rational). The 43 left are size and bound estimates. `Compile[{{x, _Real}}, x + 1][(10^400 + 1)/(3*10^399)]` was left unevaluated; it is `4.333333333333334` |
| F2 | about 87 `mpz_fits_slong_p` / `mpz_get_si` sites test a 32-bit `long` on LLP64 Windows | fixed (phase 7): demotions to Integer and int64 readers use `expr_mpz_get_int64`; helpers returning `long` decline a value the cast would truncate (`Zeta[2^32 + 2]` read as `Zeta[2]`). The opposite direction was worse: `expr_to_mpz`, `expr_mpz_borrow` and `expr_new_bigint_from_int64` used `mpz_init_set_si`, so every Integer above 2^31 converted wrong on Windows; they and 96 other int64 sites now use `mth_mpz_set_i64`. Open: the Integers `Solve` guards stay 32-bit until the solver's own `_si` checking arithmetic (`si_verify`) is 64-bit, or widening them would certify false solutions; about 250 arithmetic `_si`/`_ui` calls are unclassified |
| F3 | `cpow` cleanup drops an imaginary part below 1e-14 (`(1. + 1.*^-15 I)^0.5`) | kept: within the machine error |
| F4 | `checked_int.h` cites `tests/test_compile_overflow.c`, which does not exist | fixed (phase 7): the test exists and is built twice (`compile_overflow_tests` with the builtins, `compile_overflow_fallback_tests` with `-DCOMPILE_NO_OVERFLOW_BUILTIN`); each checks the overflow flag and the wrapped value of every helper, in both families, against GMP over the int32/int64 boundary set |
| F5 | no automated check guards the machine boundary across surfaces | open: an `ND`/Compile boundary sweep (see Phase 8 in `tasks/todo.md`) |
| F6 | real `Binomial`: `Binomial[1.*^8, 2.]`, `Binomial[-5., 3.]`, `Binomial[5., 7.]` were `Indeterminate` (a Gamma factor past MPFR's range; Gamma poles whose ratio is finite) | fixed (phase 7): a log-Gamma quotient on overflow; integer-valued operands use the exact path's reflection/symmetry/zero conventions, at machine and arbitrary precision |

## Known, not changed
- A promoted value does not demote back when a later result fits again:
  `(1.*^308 * 10.)/1.*^10` stays a 53-bit MPFR `1.*10^299`.
- `(1. + 1.*^-15 I)^0.5` drops its 5*10^-16 imaginary part (cpow cleanup).
  That part is within the machine error of `cpow`.
