# Integrate

!!! warning "Status: Partial"
    implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## Description

```text
Integrate[f, x] gives the indefinite integral of f with respect to x.
Integrate[f, x, Method -> "<name>"] dispatches directly to a single
subroutine, bypassing the default cascade.  Accepted method names:
  "Automatic"          — try BronsteinRational, then RischNorman, then CRCTable (default)
  "BronsteinRational"  — Integrate`BronsteinRational (polynomial / rational)
  "DerivativeDivides"  — Integrate`DerivativeDivides (substitution u(x); direct + Eliminate/Solve)
  "LinearRadicals"     — Integrate`LinearRadicals (rationalise radicals of a x + b)
  "QuadraticRadicals"  — Integrate`QuadraticRadicals (Euler substitution for Sqrt[a x^2 + b x + c])
  "LinearRatioRadicals" — Integrate`LinearRatioRadicals (rationalise radicals of (a x + b)/(c x + d))
  "Weierstrass"        — Integrate`Weierstrass (continuous tan(x/2) / tanh(x/2) substitution)
  "RischNorman"        — Integrate`RischNorman (Bronstein pmint heuristic)
  "CRCTable"           — Integrate`CRCTable (lazy-loaded CRC integral table)
  "Undefined"          — Integrate`Undefined (unknown functions u[x], u'[x]; Roach §1.7)
Named methods are strict: failure returns unevaluated, with no fallback.
The CRCTable rules are loaded from disk on first use only.
An applied 1-D InterpolatingFunction integrates to its antiderivative
InterpolatingFunction (mirroring D).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Integrate[3 + 5 x + 2 x^2, x]
Out[1]= 3 x + 5/2 x^2 + 2/3 x^3

In[2]:= Integrate[2 x/(x^2 + 1), x]
Out[2]= Log[1 + x^2]

In[3]:= Integrate[1/(x - a)^2, x]
Out[3]= -1/(-a + x)

In[4]:= Integrate[(2x+3)/(x^2+3x+5)^2, x]
Out[4]= -1/(5 + 3 x + x^2)

In[5]:= Integrate[1/((x-1)(x-2)(x-3)), x]          (* Phase 2 LRT closes this *)
Out[5]= -Log[-2 + x] + 1/2 Log[3 - 4 x + x^2]

In[6]:= Integrate[1/(x^2 + 1), x]                  (* Phase 4 LogToReal *)
Out[6]= ArcTan[x]

In[7]:= Integrate[1/(x^4 + x^2 + 1), x]            (* two quadratic factors *)
Out[7]= 1/4 Log[1 + x + x^2] + 1/2 ArcTan[(-1 + 2 x)/Sqrt[3]]/Sqrt[3] + 1/2 ArcTan[(1 + 2 x)/Sqrt[3]]/Sqrt[3] - 1/4 Log[1 - x + x^2]

In[8]:= Integrate[Sin[x], x, Method -> "RischNorman"]  (* strict, no fallback *)
Out[8]= -Cos[x]
```

## Implementation notes

**Algorithm.** `Integrate[f, x]` computes an *indefinite, single-variable*
antiderivative only; there is no constant of integration and no definite-bound
support, and results are not guaranteed to be simplified. The dispatcher
`builtin_integrate` (src/calculus/integrate.c) routes through a fixed cascade of
context-qualified sub-integrators, each of which returns the antiderivative on
success or comes back as its own unevaluated head (e.g.
`Integrate\`BronsteinRational[...]`) to signal "fall through". The Automatic
cascade order is:

1. **Undefined-function integrator** (`try_undefined` → src/calculus/integrate_unknown.c),
   Roach 1992 §1.7: integrands rational in unknown `u[x]` and their derivatives
   (e.g. `x f'[x] + f[x] -> x f[x]`). Cheaply gated.
2. **Bronstein rational integration** (`try_rational` → src/calculus/intrat.c),
   gated by `PolynomialQ`/rational-function tests. This is the mathematically
   substantial stage: polynomial part by term-wise power rule, then **Mack's
   linear Hermite reduction** (`HermiteReduce`) to split off the rational
   (algebraic) part of the antiderivative, then the **Lazard–Rioboo–Trager
   log part** (`IntRationalLogPart`) built on a subresultant PRS, with Rioboo's
   recursive `LogToAtan`/`LogToReal` conversions producing real-elementary
   `Log`/`ArcTan` forms (a `RootSum`-based `NaiveLogPart` is the universal
   fallback). A pre-Hermite derivative-recognition fast path catches `c·D'/D^k`.
3. **Radical substitutions**: linear-radical (`try_linrad`), quadratic-radical
   Euler substitution (`try_quadrad`), and linear-ratio Möbius substitution
   (`try_linratiorad`).
4. **Derivative-divides** substitution (`try_derivdivides`): recognises
   `c·h(u(x))·u'(x)` and reduces to `Integrate[h[u], u]`.
5. **Risch–Norman heuristic** (`try_risch` → src/calculus/intrischnorman.c),
   Bronstein's *pmint*: the parallel-Risch ansatz for transcendental
   integrands (exp/log/trig via `tan` rewriting). Builds a candidate monomial
   set, a vector field with `splitFactor`/deflation, a linear system solved by
   `RowReduce`, plus a log-candidate sum; retries over `K = I`.
6. **CRC integral table** (`try_crctable`): a large rule set lazily
   `Get`-loaded from src/internal/CRCMathTablesIntegrals.m on first use only.
   Recursion-capped at `MAX_CRC_DEPTH` (256); a leaked internal
   `IntegrateTable[...]` head anywhere in the result counts as a miss.

An explicit `Method -> "<name>"` option (parsed by `parse_method_option`)
bypasses the cascade and dispatches to one stage strictly (no fallback).
Inexact integrands are force-rationalised by `common_rationalize_input`
(min-precision-aware, ½-ulp fallback), integrated exactly, then numericalised
back via `common_numericalize_result` for inexact-in/inexact-out semantics.
Applied 1-D `InterpolatingFunction` objects integrate to a fresh antiderivative
`InterpolatingFunction` before any rationalisation (`integrate_interp`).

**Data structures.** Everything is a Mathilda `Expr*` tree; the rational stage
leans on Mathilda's own polynomial builtins (`Together`, `Numerator`,
`Denominator`, `SubresultantPolynomialRemainders`, `PolynomialGCD`). Stages
communicate by evaluating context-qualified heads (`call_stage`) and
distinguishing success from passthrough by inspecting the result head
(`result_is_unresolved`, `result_contains_head`).

**Complexity / limits.** Indefinite, single-variable only. The rational stage
is complete for rational integrands (Hermite + LRT); transcendental coverage is
heuristic (pmint may give up) backed by a finite table. No definite integrals,
no multivariate integration, no constant of integration.

- `Protected`, `Listable`.
- Nine-stage dispatch cascade (`DerivativeDivides`, `LinearRadicals`,

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Partial** — implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## References

- Bronstein, "Symbolic Integration I: Transcendental Functions", 2nd ed. (Springer, 2005).
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (Kluwer, 1992), ch. 11–12.
- M. Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (Springer, 2005).
- M. Bronstein, "poor man's integrator" (pmint), 2004 — parallel Risch / Risch–Norman heuristic.
- B. M. Trager, "Algebraic Factoring and Rational Function Integration", ACM SYMSAC 1976.
- K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992).
- K. Roach, "Symbolic Integration", 1992 (undefined-function integrands, §1.7).
- CRC Standard Mathematical Tables, 31st ed. — reference integral tables.
- Source: [`src/calculus/integrate.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/integrate.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Integrate[1/(1 + x^2), x]
Out[1]= ArcTan[x]
```

```mathematica
In[1]:= Integrate[1/x, x]
Out[1]= Log[x]
```

```mathematica
In[1]:= Integrate[Cos[x], x]
Out[1]= Sin[x]
```

```mathematica
In[1]:= Integrate[x^3 + x, x]
Out[1]= 1/2 x^2 + 1/4 x^4
```

### Notes

`Integrate[f, x]` computes the indefinite integral via a cascade: Bronstein's rational-function algorithm, then the Risch–Norman (`pmint`) heuristic, then the lazy-loaded CRC integral tables; `Method -> "<name>"` pins a single subroutine. Antiderivatives are returned without an integration constant and are not always simplified — for example `Integrate[Sin[x], x]` returns `-(1 + Cos[x])` rather than `-Cos[x]`. Definite integration is **not** supported in this build: `Integrate[x^2, {x, 0, 1}]` threads `Integrate` over the bound list and returns the garbage form `{1/3 x^3, Integrate[x^2, 0], Integrate[x^2, 1]}` instead of `1/3`. Restrict use to the indefinite, single-variable form.
