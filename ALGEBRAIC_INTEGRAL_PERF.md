---
title: "Algebraic integration stress test: Trager and ParallelMixedTower on *Diabolically Difficult Integrals*, Vol. 1"
date: 2026-10-30 (Mathilda v0.6087)
---

# Summary

Two of Mathilda's algebraic integrators were run over every integral of S. Blake, *Diabolically Difficult Integrals*,
Volume 1, with FriCAS as a third opinion:

- `Integrate[f, x, Method -> "Trager"]`: Trager's algorithm for a rational function of x and one radical over Q.
- `Integrate[f, x, Method -> "ParallelMixedTower"]` (PMT): the parallel Risch method over mixed towers.

There were 803 integrands. 207 have symbolic parameters, and each of those also appears with its parameters set to
small integers.

**The main results**

1. **Trager is sound throughout.** It gave no wrong answers and no false non-elementarity certificates. It finds 338
   antiderivatives and 15 proofs of non-elementarity. In 113 cases it gets an answer where FriCAS does not.
   - Its median time is 0.64 s.
   - Its limits come from its scope. It declines nested radicals (207), symbolic parameters (109 at the entry gate),
     several radicals (51) and transcendental integrands (52).
2. **PMT covers more classes but had soundness defects.**
   - It returned a leaked internal symbol as the "antiderivative" in 56 cases. This is fixed.
   - It issued 2 false certificates, through a precision bug in `Series`. A fix is in hand.
   - It returned 13 answers that are right only on part of the domain.
   - One case crashed the kernel through a FLINT exception.
   - 201 cases time out at 120 s.
3. **Together the two methods decide 476 integrands.** FriCAS decides 528 at a 120 s timeout.
   - 130 integrands are decided by Mathilda but not FriCAS.
   - 182 are decided by FriCAS but by neither Mathilda method.
   - 145 are decided by no system.
4. **On your parameter question.** Trager decides 86 of the 207 parametric integrals once the parameters are
   instantiated: 83 antiderivatives and 3 proofs. Extending Trager to Q(a, b, …) should win most of those 86.
   PMT already accepts parameters but decides only 29 symbolic cases. FriCAS solves 123.

# Method

**Corpus**

- File: `tests/data/ddi_vol1_corpus.m`, one record per line: id, chapter, integrand, expected status, book reference.
- Transcribed from the LaTeX chapter files (`*_problems.tex`, with `*_hints.tex` and `*_solutions.tex` used for the
  expected status).
- Families over n are expanded at the listed n.
- Skipped: definite integrals, non-integration exercises, ArcSin/ArcCsc integrands, and templates without a concrete
  integrand.
- **Audit:** every record was checked against the LaTeX, symbol by symbol. No integrand errors were found. 16 expected
  labels were corrected to say only what the book states.

**Harness**

- Tools: `tools/ddi_stress.py` runs the cases and `tools/ddi_report.py` writes this report.
- Each (integrand, method) pair runs in a fresh process under `TimeConstrained[…, 120]`.
- With an explicit `Method -> "ParallelMixedTower"`, PMT runs with no internal time budget (changed in v0.6087).
- A Mathilda antiderivative is checked two ways:
  - `PossibleZeroQ[D[F, x] - f]`;
  - a 40-digit numeric check at three real points where the integrand is real, with the parameters instantiated.
- FriCAS answers are checked numerically.
- A FriCAS answer in special functions (erf, dilog, …) counts as FriCAS deciding the integral is not elementary.

**Status labels**

- WRONG: the derivative check fails.
- BRANCH: D[F] = f on part of the real line and −f elsewhere.
- DISAGREE: a non-elementarity certificate for an integral that FriCAS or the book shows to be elementary.
- BOOK-CONFLICT: a verified antiderivative for an integral the book calls non-elementary.

**Caveat about verification.** These baseline verdicts accept an answer that passes the numeric check even when
`PossibleZeroQ` says False. That happened for 24 Trager and 54 PMT answers. The bar is now higher on two counts:

- An antiderivative must be valid for all complex x, not just on the real line.
- `PossibleZeroQ` alone decides. A `PossibleZeroQ` false negative on a correct answer is a `PossibleZeroQ` bug to fix,
  not something to patch over with numerics.

Under that bar some of these 78 answers will turn out to be branch errors, and the rest are `PossibleZeroQ` false
negatives. Both are in the fix plan below.

# Defects found

| # | defect | component | status |
|-|----------------------------------------|--------------|---------|
| 1 | PMT returned the unassigned local ``ParallelMixed`Private`r$NNN`` as an antiderivative (56 cases, also through the default `Integrate`). `Return[False]` inside a `With` in `eulerStep` left only the `With`; a self-referential `Dwp` then overflowed `$RecursionLimit`. | `ParallelMixed.m` (BuildTower) | **fixed** (v0.6087), regression test |
| 2 | False certificates: `Sqrt[1-x] Sqrt[1+x] Sqrt[1-x^4]` and multiple_radicals_4 (iii) "proved" non-elementary. `Series` came back short of the requested order at poles of order ≥ 7 (fixed padding), so PMT read a residue of −5/4 as 0. | `src/calculus/series.c` | fix in hand (Phase 1a) |
| 3 | The verify gate let a recursion-overflow result through. After `$RecursionLimit` the evaluator runs degraded, and the gate's decline branches were skipped. `fac` is a global variable. | `ParallelMixed.m`, `integrate.c` | planned (Phase 1b) |
| 4 | A kernel abort (FLINT `qqbar_pow_fmpq: division by zero`) on `((2x^3-10) Sqrt[2x^3+3x^2+5])/(2x^3+5)^2` | `src/poly/flint_qqbar.c` | planned (Phase 1c) |
| 5 | 13 answers right on one real interval and of the wrong sign on another: separate radicals are rebuilt as one principal radical, and a square factor's sign is fixed at one sample point | `ParallelMixed.m` (`rootNormalize`, `signFix`, `branchUnit`) | planned (Phase 2b) |
| 6 | 78 verified answers where `PossibleZeroQ` says False | PossibleZeroQ / both methods | planned (Phase 2d) |
| 7 | PMT: 201 timeouts and 151 "tower construction failed" declines (one radical only, no u = x^k substitution, conics with no rational point) | `ParallelMixed.m` | planned (Phase 3) |
| 8 | Trager: rational functions (11) and algebraic constants (8) declined; no parameters (207) | `src/internal/trager/` | planned (Phase 4) |
| 9 | Trager: nested radicals (207) and several radicals (51) | `src/internal/trager/` | separate plan (general function fields) |

**Harness defects found and fixed along the way:**

- An unverifiable parametric answer was counted as WRONG. Parameters are now instantiated before differentiating.
- FriCAS answers in special functions were counted as elementary.
- The PMT decline reasons were parsed wrongly.

# The parameter question

- AXIOM and FriCAS integrate over constant fields Q(a, b, …). For the torsion and reduction steps they specialise the
  parameters.
- Mathilda's Trager works over Q only. That is not just the entry gate: the constant fields in `Field.m`, the residue
  fields in `Residues.m`, and the reductions mod p in `TorsionModP.m`, `Genus.m` and `Jacobian.m` all assume Q.
- On this corpus, 86 of the 207 parametric integrals are decided by Trager once instantiated. So extending Trager to
  Q(params) has a measured payoff, and it is Phase 4b of the plan.
- One soundness rule carries over. A non-elementarity verdict must be proved generically; specialisation may only
  supply bounds and orders.

# Fix plan (in order)

1. **Soundness.**
   - The `Series` padding bug.
   - The PMT gate accepts only on positive evidence.
   - FLINT exceptions never abort the kernel.
2. **Answers valid for all complex x.**
   - Complex sample points in the harness.
   - PMT back-substitutes in the integrand's own radicals.
   - Power-free radicands.
   - `PossibleZeroQ` becomes the sole verifier, with its false negatives fixed.
3. **PMT coverage and speed.**
   - A u = x^k substitution pre-pass.
   - Euler steps through algebraic points.
   - Partially ramified places.
   - Profile and fix the timeouts by algorithm, with no time caps.
4. **Trager.**
   - Route rational integrands and algebraic constants.
   - Work over Q(params).
   - Fix the 4 timeouts.
5. **Separate plan:** Trager for general algebraic function fields (primitive element and integral basis), for nested
   and multiple radicals.

# Findings about the book

- ∫ ArcTan[x²] dx is elementary: `x ArcTan[x^2] - Integrate[2 x^2/(1 + x^4), x]`. The book (inverse_function_3) says it is not.
  Mathilda's verified antiderivative is listed as a book conflict.
- Possible typos:
  - realis_2 (v): the sign pairing contradicts special case (i);
  - realis_2 (iv);
  - multiple_radicals_0 (vi) has a stray constant `+(a+b)`;
  - timofeev_3 (iv): a stale LaTeX comment says 1/3 where the text says 3/2;
  - multiple_methods problem 1: the LaTeX comment differs from the displayed integral.
- Some expected labels rest on the chapter's framing rather than an explicit statement: generalised_gunther and
  multiple_radicals have empty hint and solution files. The book's own "elementary" counts above use these labels.

# Results

<!-- BEGIN GENERATED -->
Run of 2026-10-30 15:43:22: **803 integrands** (207 with symbolic parameters, 596 over Q), per-method timeout 120 s, FriCAS timeout 120 s, wall time 76 min; ParallelMixedTower re-run on 803 cases (33 min).

### Headline

|  | Trager | ParallelMixedTower | FriCAS |
|---|---|---|---|
| verified antiderivative | 338 | 298 | 482 |
| proved non-elementary | 15 | 22 | 46 |
| declined | 446 | 268 | 0 |
| timeout | 4 | 201 | 253 |
| error / crash | 0 | 1 | 21 |
| wrong answer | 0 | 8 | 1 |
| right on part of the real line only (branch) | 0 | 5 | 0 |
| answer not verifiable in budget | 0 | 0 | 0 |

Of the **646** integrands the book states are elementary: Trager solves **281**, ParallelMixedTower solves **258**, FriCAS **419**, either Mathilda method **381**.

- over Q (no parameters) (596): Trager solved 338, decided 353; ParallelMixedTower solved 269, decided 291; FriCAS solved 359.
- with symbolic parameters (207): Trager solved 0, decided 0; ParallelMixedTower solved 29, decided 29; FriCAS solved 123.

### Classification against FriCAS and the book

| class | Trager | ParallelMixedTower |
|---|---:|---:|
| agree elem | 225 | 257 |
| Mathilda only | 113 | 40 |
| agree non-elem | 9 | 14 |
| non-elem unconf. | 6 | 6 |
| FriCAS only | 257 | 210 |
| neither | 190 | 142 |
| timeout | 3 | 117 |
| error | 0 | 1 |
| WRONG | 0 | 8 |
| wrong branch on part of R | 0 | 5 |
| DISAGREE | 0 | 2 |
| book conflict | 0 | 1 |

### Per chapter

Cells are *verified antiderivatives / proved non-elementary / declined / timeout+error+wrong*.

| chapter | n | book elem | Trager | ParallelMixedTower | FriCAS (elem/non) |
|---|---:|---:|---|---|---|
| warmup | 73 | 73 | 35 / 0 / 37 / 1 | 53 / 0 / 9 / 11 | 73 / 0 |
| laurent_poly_subst | 60 | 60 | 43 / 0 / 17 / 0 | 20 / 0 / 9 / 31 | 42 / 0 |
| rational_subs | 32 | 32 | 27 / 0 / 5 / 0 | 15 / 0 / 3 / 14 | 16 / 0 |
| gunther | 33 | 29 | 21 / 4 / 8 / 0 | 14 / 4 / 7 / 8 | 26 / 0 |
| power_expand | 6 | 6 | 5 / 0 / 1 / 0 | 5 / 0 / 0 / 1 | 6 / 0 |
| generalised_gunther | 88 | 80 | 66 / 2 / 20 / 0 | 19 / 0 / 11 / 58 | 19 / 2 |
| chebychev | 87 | 51 | 3 / 7 / 77 / 0 | 17 / 7 / 47 / 16 | 32 / 32 |
| inverse_function | 18 | 10 | 0 / 0 / 18 / 0 | 2 / 5 / 10 / 1 | 6 / 5 |
| product_rule | 21 | 21 | 19 / 0 / 2 / 0 | 21 / 0 / 0 / 0 | 21 / 0 |
| log_term | 16 | 16 | 16 / 0 / 0 / 0 | 16 / 0 / 0 / 0 | 16 / 0 |
| realis | 22 | 7 | 17 / 2 / 3 / 0 | 11 / 2 / 8 / 1 | 18 / 0 |
| nested_radicals | 161 | 159 | 0 / 0 / 161 / 0 | 40 / 1 / 98 / 22 | 86 / 2 |
| multiple_radicals | 55 | 34 | 0 / 0 / 55 / 0 | 9 / 2 / 26 / 15 | 40 / 2 |
| euler | 22 | 22 | 10 / 0 / 11 / 1 | 6 / 0 / 9 / 7 | 13 / 0 |
| timofeev | 23 | 7 | 19 / 0 / 4 / 0 | 15 / 0 / 3 / 5 | 23 / 0 |
| potpourri | 43 | 7 | 31 / 0 / 12 / 0 | 22 / 0 / 13 / 6 | 27 / 0 |
| mult_piecewise | 1 | 0 | 0 / 0 / 1 / 0 | 0 / 0 / 1 / 0 | 0 / 0 |
| danger_zone | 20 | 18 | 17 / 0 / 1 / 2 | 1 / 0 / 9 / 10 | 4 / 0 |
| multiple_methods | 4 | 0 | 2 / 0 / 2 / 0 | 0 / 0 / 1 / 3 | 0 / 0 |
| introduction | 18 | 14 | 7 / 0 / 11 / 0 | 12 / 1 / 4 / 1 | 14 / 3 |

### Timing (seconds, cases each system answered — antiderivative or proof)

| system | answered | median | mean | max |
|---|---:|---:|---:|---:|
| Trager | 353 | 0.644 | 1.774 | 31.42 |
| ParallelMixedTower | 320 | 0.273 | 3.227 | 119.55 |
| FriCAS | 528 | 0.200 | 2.452 | 90.90 |

On the 185 integrands both methods solve, median LeafCount of the antiderivative: Trager 118, ParallelMixedTower 101.

### Why each method declines

**Trager** (446 declines)

| reason | count |
|---|---:|
| nested radicals | 207 |
| radicand is not a rational function over Q | 109 |
| transcendental function in the integrand | 52 |
| more than one radicand (a radical compositum) | 51 |
| no radical in the integrand (a rational function) | 11 |
| algebraic constant in the integrand (coefficients must be rational) | 8 |
| integrand is not a rational function of x and the radical over Q | 6 |
| integrand is not algebraic over Q(x) | 2 |

**ParallelMixedTower** (268 declines)

| reason | count |
|---|---:|
| tower construction failed | 151 |
| no solution within bounds | 46 |
| needs torsion realisation (milestone iii) | 29 |
| partially ramified branch prime | 15 |
| residue class lost | 10 |
| verification failed | 6 |
| certificate withheld | 5 |
| internal | 4 |
| deep residue sheet undecided | 2 |

### Wrong answers (failed derivative check) (8)

| id | method | integrand | book ref | detail |
|-----------|------|------------------|---------|--------------|
| multrad\-intro\-1 | PMT | \(x\^4\-12\)/\(\(x\^4\-8\*x\^2\+12\)\*Sqrt\[x\^2\-4\]\*Sqrt\[x\^2\-3\]\) | Section intro example \(solved, \-ArcTanh\) | D\[F\] = \-f at every sample \(global sign\) num=\[\['7/3', 'neg', '9.812729'\], \['13/5', 'neg', '5.783089'\], \['11/3',… |
| multrad\-multiple\_radicals\_0\-i | PMT | \(x\^4\-a\*b\)/\(\(d\*x\^4\-\(a\*d\+b\*d\-c\)\*x\^2\+a\*b\*d\)\*Sqrt\[x\^2\-a\]\*Sqrt\[x\^2\-b\]\) | Problem multiple\_radicals\_0 \(i\) | D\[F\] = \-f at every sample \(global sign\) num=\[\['7/3', 'neg', '0.3868618'\], \['13/5', 'neg', '0.2779431'\], \['11/3… |
| multrad\-multiple\_radicals\_0\-i\-inst | PMT | \(x\^4\-2\*3\)/\(\(7\*x\^4\-\(2\*7\+3\*7\-5\)\*x\^2\+2\*3\*7\)\*Sqrt\[x\^2\-2\]\*Sqrt\[x\^2\-3\]\) | Problem multiple\_radicals\_0 \(i\); inst a=… | D\[F\] = \-f at every sample \(global sign\) num=\[\['7/3', 'neg', '0.1891279'\], \['13/5', 'neg', '0.1179702'\], \['11/3… |
| multrad\-multiple\_radicals\_0\-ii | PMT | \(x\^4\-a\*b\)/\(\(x\^4\+a\*b\)\*Sqrt\[x\^2\-a\]\*Sqrt\[x\^2\-b\]\) | Problem multiple\_radicals\_0 \(ii\) | D\[F\] = \-f at every sample \(global sign\) num=\[\['7/3', 'neg', '0.4375095'\], \['13/5', 'neg', '0.346314'\], \['11/3'… |
| multrad\-multiple\_radicals\_0\-ii\-inst | PMT | \(x\^4\-2\*3\)/\(\(x\^4\+2\*3\)\*Sqrt\[x\^2\-2\]\*Sqrt\[x\^2\-3\]\) | Problem multiple\_radicals\_0 \(ii\); inst a… | D\[F\] = \-f at every sample \(global sign\) num=\[\['7/3', 'neg', '0.4571958'\], \['13/5', 'neg', '0.3630166'\], \['11/3… |
| multrad\-multiple\_radicals\_0\-iv | PMT | \(x\^4\-a\*b\)\*Sqrt\[x\^2\-a\]\*Sqrt\[x\^2\-b\]/\(x\^2\*\(d\*x\^4\+\(c\-a\*d\-b\*d\)\*x\^2\+a\*b\*d\)\) | Problem multiple\_radicals\_0 \(iv\) | D\[F\] = \-f at every sample \(global sign\) num=\[\['7/3', 'neg', '1.058826'\], \['13/5', 'neg', '1.101532'\], \['11/3',… |
| multrad\-multiple\_radicals\_0\-iv\-inst | PMT | \(x\^4\-2\*3\)\*Sqrt\[x\^2\-2\]\*Sqrt\[x\^2\-3\]/\(x\^2\*\(7\*x\^4\+\(5\-2\*7\-3\*7\)\*x\^2\+2\*3\*7\)\) | Problem multiple\_radicals\_0 \(iv\); inst a… | D\[F\] = \-f at every sample \(global sign\) num=\[\['7/3', 'neg', '0.2924838'\], \['13/5', 'neg', '0.312335'\], \['11/3'… |
| multrad\-multiple\_radicals\_3\-ii | PMT | \(x\^2\-2\*a\*x\+a\*b\+a\*c\-b\*c\)/\(\(x\^2\-\(b\+c\+1\)\*x\+a\+b\*c\)\*Sqrt\[x\-a\]\*Sqrt\[x\-b\]\*Sqrt\[x\-c\]\) | Problem multiple\_radicals\_3 \(ii\) | D\[F\] = \-f at every sample \(global sign\) num=\[\['7/3', 'neg', '5.160072'\], \['13/5', 'neg', '8.15654'\], \['11/3', … |

### Right on part of the real line only (D[F] = f on one interval, -f on another) (5)

| id | method | integrand | book ref | detail |
|-----------|------|------------------|---------|--------------|
| multrad\-multiple\_radicals\_3\-i\-inst | PMT | \(x\^2\-2\*2\*x\+2\*3\+2\*5\-3\*5\)\*Sqrt\[x\-3\]\*Sqrt\[x\-5\]/\(\(x\^2\-\(3\+5\-1\)\*x\+3\*5\-2\)\*\(x\-2\)\^\(3/2\)\) | Problem multiple\_radicals\_3 \(i\); inst a=… | D\[F\] = f at 1 point\(s\), \-f at 2; points \[\['23/3', 'ok', '0.0'\], \['3/7', 'neg', '0.1813736'\], \['5/4', 'neg', '3… |
| multrad\-multiple\_radicals\_3\-ii\-inst | PMT | \(x\^2\-2\*2\*x\+2\*3\+2\*5\-3\*5\)/\(\(x\^2\-\(3\+5\+1\)\*x\+2\+3\*5\)\*Sqrt\[x\-2\]\*Sqrt\[x\-3\]\*Sqrt\[x\-5\]\) | Problem multiple\_radicals\_3 \(ii\); inst a… | D\[F\] = f at 1 point\(s\), \-f at 2; points \[\['23/3', 'ok', '0.0'\], \['3/7', 'neg', '0.01852804'\], \['5/4', 'neg', '… |
| multrad\-multiple\_radicals\_3\-iii\-inst | PMT | \(x\^2\-2\*3\*x\+2\*3\-2\*5\+3\*5\)\*Sqrt\[x\-3\]/\(\(x\^2\-\(2\+5\-1\)\*x\-3\+2\*5\)\*\(x\^2\-\(2\+5\+1\)\*x\+3\+2\*5\)\*Sqrt\[x\-2\]\*S… | Problem multiple\_radicals\_3 \(iii\); inst … | D\[F\] = f at 1 point\(s\), \-f at 2; points \[\['23/3', 'ok', '0.0'\], \['3/7', 'neg', '0.2290425'\], \['5/4', 'neg', '1… |
| potpourri\-p7\-ii | PMT | 1/\(\(x\+1\)\*\(x\^4\-x\^2\-2\*x\+2\)\^\(3/2\)\) | Potpourri problem 7 \(ii\) | D\[F\] = f at 1 point\(s\), \-f at 2; points \[\['3/7', 'ok', '0.0'\], \['5/4', 'neg', '3.811096'\], \['7/3', 'neg', '0.0… |
| potpourri\-p8\-i | PMT | Sqrt\[1\-x\^2\]/Sqrt\[1\-Sqrt\[1\-x\^2\]\] | Potpourri problem 8 \(nested radicals\) \(i… | D\[F\] = f at 2 point\(s\), \-f at 1; points \[\['3/7', 'ok', '0.0'\], \['2/9', 'ok', '0.0'\], \['\-1/3', 'neg', '7.884789… |

### Antiderivatives not verifiable within the budget (0)

None.

### Disagreements (proved non-elementary, but FriCAS or the book has an antiderivative) (2)

| id | method | integrand | book ref | detail |
|-----------|------|------------------|---------|--------------|
| multrad\-multiple\_radicals\_4\-i | PMT | Sqrt\[1\-x\]\*Sqrt\[x\+1\]\*Sqrt\[1\-x\^4\] | Problem multiple\_radicals\_4 \(i\) | \{"not elementary", "holomorphic remainder: residual second\-kind differential is not exact \(every bound in forc… |
| multrad\-multiple\_radicals\_4\-iii | PMT | Sqrt\[x\+1\]\*Sqrt\[1\-x\^4\]/Sqrt\[1\-x\] | Problem multiple\_radicals\_4 \(iii\) | \{"not elementary", "holomorphic remainder: residual second\-kind differential is not exact \(every bound in forc… |

### Book conflicts (verified antiderivative where the book says non-elementary) (1)

| id | method | integrand | book ref | detail |
|-----------|------|------------------|---------|--------------|
| inverse\_function\-inverse\_function\_3\-iii | PMT | ArcTan\[x\^2\] | Problem inverse\_function\_3 \(iii\) BOOK CL… | \(\-ArcTan\[x\^2\] \+ Sqrt\[2\] x ArcTan\[x\^2\]\)/Sqrt\[2\] \- Sqrt\[2\] ArcTan\[1 \+ Sqrt\[2\] x\] \- Log\[1 \- Sqrt\[2\] x \+ x\^2\]/\(2 S… |

### Crashes / errors (1)

| id | method | integrand | book ref | detail |
|-----------|------|------------------|---------|--------------|
| laurent\_poly\_subst\-prob\_laurent\_sub\_13\-i… | PMT | \(\(2\*x\^3\-2\*5\)\*Sqrt\[2\*x\^3\+3\*x\^2\+5\]\)/\(2\*x\^3\+5\)\^2 | Problem 13 \(ii\) instantiated a=2,b=3,c=5 | rc=\-6; FLINT exception \(Divide by zero\):     qqbar\_pow\_fmpq: division by zero |

### Timeouts (205)

| id | method | integrand | book ref | detail |
|-----------|------|------------------|---------|--------------|
| warmup\-prob\_warmup\_4\-iii | PMT | Sqrt\[a\*x\^2\+b\*x\+c\]/x\^2 | Problem 0.4 \(iii\) Euler substitution | FriCAS ELEM; book elem |
| warmup\-prob\_warmup\_4\-iv | PMT | x\^2\*Sqrt\[a\*x\^2\+b\*x\+c\] | Problem 0.4 \(iv\) Euler substitution | FriCAS ELEM; book elem |
| warmup\-prob\_warmup\_4\-v | PMT | \(a\*x\^2\+b\*x\+c\)\^\(3/2\)/x\^2 | Problem 0.4 \(v\) Euler substitution | FriCAS ELEM; book elem |
| warmup\-prob\_warmup\_4\-vi | PMT | x\*\(a\*x\^2\+b\*x\+c\)\^\(3/2\) | Problem 0.4 \(vi\) Euler substitution | FriCAS ELEM; book elem |
| warmup\-prob\_warmup\_4\-vii | PMT | \(a\*x\^2\+b\*x\+c\)\^\(5/2\)/x | Problem 0.4 \(vii\) Euler substitution | FriCAS ELEM; book elem |
| warmup\-prob\_warmup\_4\-viii | PMT | x\*\(a\*x\^2\+b\*x\+c\)\^\(5/2\) | Problem 0.4 \(viii\) Euler substitution | FriCAS ELEM; book elem |
| warmup\-prob\_warmup\_4\-ix | PMT | 1/\(\(x\+b\)\^2\*\(a\*x\^2\+b\*x\+c\)\^\(3/2\)\) | Problem 0.4 \(ix\) Euler substitution | FriCAS ELEM; book elem |
| warmup\-prob\_warmup\_4\-x | PMT | 1/\(\(x\+b\)\^2\*\(a\*x\^2\+b\*x\+c\)\^\(5/2\)\) | Problem 0.4 \(x\) Euler substitution | FriCAS ELEM; book elem |
| warmup\-prob\_warmup\_4\-xi | PMT | x/\(\(a\*x\^2\+b\*x\+c\)\^\(3/2\)\-\(a\*x\^2\+b\*x\+c\)\^\(5/2\)\) | Problem 0.4 \(xi\) Euler substitution | FriCAS ELEM; book elem |
| warmup\-prob\_warmup\_5\-inst | Trager | \(\(2\*x\+3\)\^\(2/3\)\+1\)/\(\(2\*x\+3\)\^\(1/2\)\+\(2\*x\+3\)\^\(1/4\)\) | Problem 0.5 multiple radicals of ax\+b in… | FriCAS ELEM; book elem |
| warmup\-prob\_warmup\_6 | PMT | \(\(\(x\-1\)/\(2\*x\+1\)\)\^\(1/4\)\-3\*\(\(x\-1\)/\(2\*x\+1\)\)\^\(3/4\)\)/\(\(x\-1\)\*\(x\+1\)\^2\*\(2\*x\-1\)\) | Problem 0.6 Moebius radicand | FriCAS ELEM; book elem |
| warmup\-prob\_warmup\_13\-iii | PMT | 1/\(x\*Sqrt\[a\*x\^2\+b\*x\+c\]\) | Problem 1.13 \(iii\) | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_1 | PMT | \(a\*x\^2\-b\)/\(\(a\*x\^2\+b\)\*Sqrt\[a\^2\*x\^4\+b\^2\]\) | Problem 1 | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_3 | PMT | \(x\^12\+1\)/\(\(x\^12\-1\)\*\(x\^8\+1\)\^\(1/4\)\) | Problem 3 | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_5\-i\-… | PMT | \(\(x\^5\+1\)\*\(4\*x\^5\-1\)\)/\(x\*\(x\^5\-2\*x\+1\)\*Sqrt\[x\^6\+x\]\) | Problem 5 \(i\) instantiated a=2 | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_5\-ii… | PMT | \(\(x\^5\-1\)\*\(4\*x\^5\+1\)\)/\(x\*\(x\^5\-2\*x\-1\)\*Sqrt\[x\^6\-x\]\) | Problem 5 \(ii\) instantiated a=2 | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_7\-i | PMT | \(c\*x\^3\-d\)/\(a\*x\^3\-b\)\^\(1/3\) | Problem 7 \(i\) | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_7\-i\-… | PMT | \(5\*x\^3\-7\)/\(2\*x\^3\-3\)\^\(1/3\) | Problem 7 \(i\) instantiated a=2,b=3,c=5,d… | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_7\-ii | PMT | \(c\*x\^3\+d\)/\(a\*x\^3\-b\)\^\(1/3\) | Problem 7 \(ii\) | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_7\-ii… | PMT | \(5\*x\^3\+7\)/\(2\*x\^3\-3\)\^\(1/3\) | Problem 7 \(ii\) instantiated a=2,b=3,c=5,… | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_8\-i\-… | PMT | \(2\*x\^4\+3\*x\^3\)\^\(1/4\) | Problem 8 \(i\) instantiated a=2,b=3 | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_8\-ii… | PMT | \(2\*x\^8\+3\*x\^7\)\^\(1/4\) | Problem 8 \(ii\) instantiated a=2,b=3 | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_8\-ii… | PMT | \(2\*x\^12\+3\*x\^11\)\^\(1/4\) | Problem 8 \(iii\) instantiated a=2,b=3 | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_9 | PMT | \(a\*x\^4\-b\)/\(x\^4\*\(2\*a\*x\^4\-b\)\^\(1/4\)\) | Problem 9 \(same as Problem 14b \(i\)\) | FriCAS TIMEOUT; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_10 | PMT | \(\(x\^6\-2\)\*\(x\^6\+1\)\*\(x\^6\-x\^4\+1\)\^\(1/4\)\)/\(x\^6\*\(x\^6\-2\*x\^4\+1\)\) | Problem 10 | FriCAS TIMEOUT; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_12\-i… | PMT | \(\(x\^3\-1\)\*\(x\^6\-x\^3\+1\)\*\(x\^6\+1\)\^\(2/3\)\)/\(x\^6\*\(x\^3\+1\)\) | Problem 12 \(ii\) | FriCAS TIMEOUT; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_13\-i | PMT | \(a\*x\^2\+b\)/\(\(a\*x\^2\-b\)\*Sqrt\[a\*x\^3\+b\*x\^2\-b\*x\]\) | Problem 13 \(i\) | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_13\-v… | PMT | \(\(2\*x\^5\-3\)\*Sqrt\[x\^6\+2\*x\^4\+x\]\)/\(\(x\^5\+1\)\*\(x\^5\+x\^3\+1\)\) | Problem 13 \(vi\) | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_13\-v… | PMT | \(\(2\*c\-a\*x\^3\)\*Sqrt\[a\*x\^3\+b\*x\^2\+c\]\)/\(\(a\*x\^3\+\(b\-3\)\*x\^2\+c\)\*\(a\*x\^3\+\(b\-2\)\*x\^2\+c\)\) | Problem 13 \(vii\) | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_14\-i | PMT | \(x\^2\+4\)/\(\(x\^2\-4\)\*\(x\^3\-x\^2\)\^\(1/3\)\) | Problem 14 \(i\) | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_14\-i… | PMT | \(\(x\^2\-4\)\*\(x\^3\+x\)\^\(1/3\)\)/\(x\^4\*\(x\^2\+2\)\) | Problem 14 \(ii\) | FriCAS TIMEOUT; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_14\-i… | PMT | \(\(x\^3\-1\)\*\(x\^6\+1\)\^\(1/3\)\)/\(x\^2\*\(x\^3\+1\)\) | Problem 14 \(iii\) | FriCAS TIMEOUT; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_14\-i… | PMT | \(\(x\^2\-2\)\*\(x\^3\+x\)\^\(1/3\)\)/\(x\^2\*\(x\^4\-2\*x\^2\+4\)\) | Problem 14 \(iv\) | FriCAS ERROR; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_14\-v | PMT | \(x\^3\*\(x\^2\+3\)\)/\(\(x\^2\+1\)\*\(x\^3\+x\^2\+1\)\*\(x\^3\-x\^2\-1\)\^\(1/3\)\) | Problem 14 \(v\) | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_14b\-… | PMT | \(x\^3\+4\*b\)/\(\(x\^3\+b\)\*\(a\*x\^4\-x\^3\-b\)\^\(1/4\)\) | Problem 14b \(ii\) | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_14b\-… | PMT | \(x\^3\+4\*3\)/\(\(x\^3\+3\)\*\(2\*x\^4\-x\^3\-3\)\^\(1/4\)\) | Problem 14b \(ii\) instantiated a=2,b=3 | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_14b\-… | PMT | \(\(x\^2\-2\)\*\(x\^2\-1\)\*\(x\^4\+x\^2\-1\)\^\(1/4\)\)/\(x\^6\*\(2\*x\^4\+x\^2\-1\)\) | Problem 14b \(iii\) | FriCAS TIMEOUT; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_16\-i… | PMT | \(\(a\*x\^2\-b\)/\(x\*\(a\*x\^2\+b\)\)\)\*Sqrt\[\(a\*x\^2\-c\*x\+b\)/\(a\*x\^2\+c\*x\+b\)\] | Problem 16 \(ii\) | FriCAS ELEM; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_17\-i… | PMT | \(\(2\*a\*x\^3\+b\)\*\(a\^2\*x\^6\-2\*a\*b\*x\^3\-a\*b\*x\^2\+b\^2\)\)/\(x\*\(a\*x\^3\-b\)\*\(a\^2\*x\^6\-2\*a\*b\*x\^3\+x\^2\+b\^2\)\)\*Sq… | Problem 17 \(ii\) | FriCAS TIMEOUT; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_18\-i | PMT | \(\(2\*x\^3\+1\)/\(x\*\(x\^3\+x\-1\)\)\)\*\(\(x\^3\-x\-1\)/\(x\^3\+x\-1\)\)\^\(1/3\) | Problem 18 \(i\) | FriCAS ERROR; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_18\-i… | PMT | \(\(2\*a\*x\^3\+b\)/\(x\*\(a\*x\^3\+x\-b\)\)\)\*\(\(a\*x\^3\-a\*x\-b\)/\(a\*x\^3\+b\*x\-b\)\)\^\(1/3\) | Problem 18 \(ii\) | FriCAS TIMEOUT; book elem |
| laurent\_poly\_subst\-prob\_laurent\_sub\_18\-i… | PMT | \(\(2\*2\*x\^3\+3\)/\(x\*\(2\*x\^3\+x\-3\)\)\)\*\(\(2\*x\^3\-2\*x\-3\)/\(2\*x\^3\+3\*x\-3\)\)\^\(1/3\) | Problem 18 \(ii\) instantiated a=2,b=3 | FriCAS TIMEOUT; book elem |
| rational\_subs\-prob\_rational\_sub\_1\-i | PMT | 1/\(\(x\+1\)\*\(3\*x\^2\+1\)\^\(1/3\)\) | Problem 1 \(i\) | FriCAS TIMEOUT; book elem |
| rational\_subs\-prob\_rational\_sub\_1\-ii | PMT | 1/\(\(x\-1\)\*\(3\*x\^2\+1\)\^\(1/3\)\) | Problem 1 \(ii\) | FriCAS TIMEOUT; book elem |
| rational\_subs\-prob\_rational\_sub\_1\-iii | PMT | 1/\(\(x\+1\)\*\(x\^3\+6\*x\^2\+3\*x\+2\)\^\(1/3\)\) | Problem 1 \(iii\) | FriCAS TIMEOUT; book elem |
| rational\_subs\-prob\_rational\_sub\_1\-iv | PMT | 1/\(\(x\-1\)\*\(x\^3\-6\*x\^2\+3\*x\-2\)\^\(1/3\)\) | Problem 1 \(iv\) | FriCAS TIMEOUT; book elem |
| rational\_subs\-prob\_rational\_sub\_1\-v | PMT | 1/\(\(x\+1\)\*\(\(a\+1\)\*x\^3\+3\*\(1\-a\)\*x\^2\+3\*\(a\+1\)\*x\-a\+1\)\^\(1/3\)\) | Problem 1 \(v\) | FriCAS TIMEOUT; book elem |
| rational\_subs\-prob\_rational\_sub\_1\-v\-inst | PMT | 1/\(\(x\+1\)\*\(\(2\+1\)\*x\^3\+3\*\(1\-2\)\*x\^2\+3\*\(2\+1\)\*x\-2\+1\)\^\(1/3\)\) | Problem 1 \(v\) instantiated a=2 | FriCAS TIMEOUT; book elem |
| rational\_subs\-prob\_rational\_sub\_1\-vi | PMT | 1/\(\(x\-1\)\*\(\(a\-1\)\*x\^3\+3\*\(a\+1\)\*x\^2\+3\*\(a\-1\)\*x\+a\+1\)\^\(1/3\)\) | Problem 1 \(vi\) | FriCAS TIMEOUT; book elem |
| rational\_subs\-prob\_rational\_sub\_1\-vi\-ins… | PMT | 1/\(\(x\-1\)\*\(\(2\-1\)\*x\^3\+3\*\(2\+1\)\*x\^2\+3\*\(2\-1\)\*x\+2\+1\)\^\(1/3\)\) | Problem 1 \(vi\) instantiated a=2 | FriCAS TIMEOUT; book elem |
| rational\_subs\-prob\_rational\_sub\_2\-iii | PMT | \(b\*x\^2\-2\*a\*x\-a\*b\)/\(\(x\^2\+a\)\*Sqrt\[2\*x\^4\+2\*b\*x\^3\+\(2\*a\+b\^2\)\*x\^2\+a\^2\]\) | Problem 2 \(iii\) | FriCAS ELEM; book elem |
| rational\_subs\-prob\_rational\_sub\_2\-iv | PMT | \(b\*x\^2\-2\*a\*x\-a\*b\)/\(x\*\(x\+b\)\*Sqrt\[2\*x\^4\+2\*b\*x\^3\+\(2\*a\+b\^2\)\*x\^2\+a\^2\]\) | Problem 2 \(iv\) | FriCAS ELEM; book elem |
| rational\_subs\-prob\_rational\_sub\_6\-i\-plus | PMT | 1/\(\(x\+1\)\*\(x\^3\+x\)\^\(1/4\)\) | Problem 6 \(i\) x\+1 | FriCAS TIMEOUT; book elem |
| rational\_subs\-prob\_rational\_sub\_6\-ii\-plu… | PMT | \(x\+2\)/\(\(x\^2\-1\)\*\(x\^3\+x\)\^\(1/4\)\) | Problem 6 \(ii\) x\+2 | FriCAS TIMEOUT; book elem |
| rational\_subs\-prob\_rational\_sub\_6\-i\-minu… | PMT | 1/\(\(x\-1\)\*\(x\^3\+x\)\^\(1/4\)\) | Problem 6 \(i\) x\-1 | FriCAS TIMEOUT; book elem |
| rational\_subs\-prob\_rational\_sub\_6\-ii\-min… | PMT | \(x\-2\)/\(\(x\^2\-1\)\*\(x\^3\+x\)\^\(1/4\)\) | Problem 6 \(ii\) x\-2 | FriCAS TIMEOUT; book elem |
| gunther\-prob\_gunther\_2\-ii | PMT | \(a\*x\-b\)/\(\(a\*x\+b\)\*Sqrt\[a\^2\*x\^3\+b\^2\*x\]\) | Problem 2 \(ii\) | FriCAS ELEM; book elem |
| gunther\-prob\_gunther\_3\-iv | PMT | \(a\*x\^2\+b\)/\(\(a\*x\^2\-b\)\*Sqrt\[a\*x\^3\+c\*x\^2\-b\*x\]\) | Problem 3 \(iv\) | FriCAS ELEM; book elem |
| gunther\-prob\_gunther\_3\-vi | PMT | \(x\^2\-1\)/\(\(x\^2\+1\)\*Sqrt\[a\*x\^4\+b\*x\^3\+c\*x\^2\+b\*x\+a\]\) | Problem 3 \(vi\) | FriCAS ELEM; book elem |
| gunther\-prob\_gunther\_4 | PMT | \(2\*a\*x\^4\+4\*a\*x\^3\-b\*x\+b\-2\*c\)/\(\(a\*x\^4\+x\^2\+\(b\+2\)\*x\+c\+1\)\*Sqrt\[a\*x\^4\+b\*x\+c\]\) | Problem 4 | FriCAS TIMEOUT; book elem |
| gunther\-prob\_gunther\_4\-inst | PMT | \(2\*2\*x\^4\+4\*2\*x\^3\-3\*x\+3\-2\*5\)/\(\(2\*x\^4\+x\^2\+\(3\+2\)\*x\+5\+1\)\*Sqrt\[2\*x\^4\+3\*x\+5\]\) | Problem 4 instantiated a=2,b=3,c=5 | FriCAS ELEM; book elem |
| gunther\-prob\_gunther\_5\-i | PMT | x/\(\(x\^3\+8\)\*Sqrt\[x\^3\-1\]\) | Problem 5 \(i\) | FriCAS ELEM; book elem |
| gunther\-prob\_gunther\_5\-ii\-pp | PMT | x/\(\(x\^3\+4\)\*Sqrt\[x\^3\+1\]\) | Problem 5 \(ii\) upper signs | FriCAS ELEM; book elem |
| gunther\-prob\_gunther\_5\-ii\-mm | PMT | x/\(\(x\^3\-4\)\*Sqrt\[x\^3\-1\]\) | Problem 5 \(ii\) lower signs | FriCAS ELEM; book elem |
| power\_expand\-power\_expand\_2\-ii | PMT | 1/\(x\^20\-3\*x\^16\+2\*x\^12\+2\*x\^8\-3\*x\^4\+1\)\^\(1/4\) | Problem 2 \(ii\) | FriCAS ELEM; book elem |
| generalised\_gunther\-intro\-1 | PMT | \(2\*x\^2\-1\)/\(\(x\^2\-1\)\*\(3\*x\^2\-1\)\*\(2\*x\^3\-x\)\^\(1/3\)\) | Chapter intro worked example \(solved\) | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_1\-main | PMT | \(2\*x\^2\+1\)/\(\(x\^2\+1\)\*\(3\*x\^2\+1\)\*\(2\*x\^3\-x\)\^\(1/3\)\) | Problem gg1 \(unlabelled, 1st of chapter\) | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_2\-i | PMT | 1/\(x\*\(x\^2\-3\*x\+2\)\^\(1/3\)\) | Problem gg2 \(i\) | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_2\-ii | PMT | \(x\-2\*a\+b\)/\(\(x\^2\+\(a\*b\^2\-2\*b\)\*x\+b\^2\-a\^2\*b\^2\)\*\(\(x\-a\)\*\(x\-b\)\)\^\(1/3\)\) | Problem gg2 \(ii\) | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_2\-ii\-inst | PMT | \(x\-2\*2\+3\)/\(\(x\^2\+\(2\*3\^2\-2\*3\)\*x\+3\^2\-2\^2\*3\^2\)\*\(\(x\-2\)\*\(x\-3\)\)\^\(1/3\)\) | Problem gg2 \(ii\) instantiated a=2,b=3 | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_4\-main | PMT | \(x\+1\)/\(x\*\(x\^4\-4\*x\^3\+7\*x\^2\-4\*x\+1\)\^\(1/4\)\) | Problem gg4 | FriCAS ELEM; book elem |
| generalised\_gunther\-prob\_gg\_5\-i | PMT | \(3\*x\^2\+1\)/\(\(x\^2\-1\)\*\(x\^4\+6\*x\^2\+1\)\^\(3/4\)\) | Problem gg5 \(i\) | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_5\-ii | PMT | \(x\-1\)\^2/\(\(x\+1\)\*\(x\^4\+6\*x\^2\+1\)\^\(3/4\)\) | Problem gg5 \(ii\) | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_5\-iii | PMT | \(x\+1\)\^2/\(\(x\-1\)\*\(x\^4\+6\*x\^2\+1\)\^\(3/4\)\) | Problem gg5 \(iii\) | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_6\-i | PMT | \(x\+a\)\^2/\(\(x\+b\)\*\(2\*x\^4\+4\*\(a\+b\)\*x\^3\+6\*\(a\^2\+b\^2\)\*x\^2\+4\*\(a\^3\+b\^3\)\*x\+a\^4\+b\^4\)\^\(3/4\)\) | Problem gg6 \(i\) | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_6\-i\-inst | PMT | \(x\+2\)\^2/\(\(x\+3\)\*\(2\*x\^4\+4\*\(2\+3\)\*x\^3\+6\*\(2\^2\+3\^2\)\*x\^2\+4\*\(2\^3\+3\^3\)\*x\+2\^4\+3\^4\)\^\(3/4\)\) | Problem gg6 \(i\) instantiated a=2,b=3 | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_6\-ii | PMT | \(x\-a\)\^2/\(\(x\-b\)\*\(2\*x\^4\-4\*\(a\+b\)\*x\^3\+6\*\(a\^2\+b\^2\)\*x\^2\-4\*\(a\^3\+b\^3\)\*x\+a\^4\+b\^4\)\^\(3/4\)\) | Problem gg6 \(ii\) | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_6\-ii\-inst | PMT | \(x\-2\)\^2/\(\(x\-3\)\*\(2\*x\^4\-4\*\(2\+3\)\*x\^3\+6\*\(2\^2\+3\^2\)\*x\^2\-4\*\(2\^3\+3\^3\)\*x\+2\^4\+3\^4\)\^\(3/4\)\) | Problem gg6 \(ii\) instantiated a=2,b=3 | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_7\-ix | PMT | \(x\^2\+x\+2\)\*\(x\^3\-1\)\^\(1/3\)/\(x\*\(x\^2\-1\)\^2\*\(x\^3\+x\^2\-2\*x\-3\)\) | Problem gg7 \(ix\) cube roots of cubics | FriCAS ERROR; book elem |
| generalised\_gunther\-prob\_gg\_8\-ii | PMT | \(x\-3\)/\(\(x\+2\)\*\(x\^2\+1\)\*\(x\^2\-3\)\^\(1/4\)\) | Problem gg8 \(ii\) fourth roots of quadrat… | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_8\-iv | PMT | \(x\+3\)\*\(x\^2\-3\)\^\(1/4\)/\(\(x\-2\)\*\(x\-1\)\^2\*\(x\^2\+1\)\) | Problem gg8 \(iv\) fourth roots of quadrat… | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_8\-v | PMT | \(x\-3\)\*\(x\^2\-3\)\^\(3/4\)/\(\(x\+1\)\^4\*\(x\+2\)\*\(x\^2\+1\)\) | Problem gg8 \(v\) fourth roots of quadrati… | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_8\-vii | PMT | x\^2\*\(x\+6\)/\(\(x\^3\+x\+4\)\*\(x\^2\+4\*x\)\^\(3/4\)\) | Problem gg8 \(vii\) fourth roots of quadra… | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_8\-viii | PMT | \(x\+3\)\*\(2\*x\+9\)/\(x\*\(x\^3\+x\+3\)\*\(x\^2\+3\*x\)\^\(3/4\)\) | Problem gg8 \(viii\) fourth roots of quadr… | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_8\-ix | PMT | \(x\+2\)\*\(x\+3\)/\(x\^3\*\(x\^3\+x\+2\)\*\(x\^2\+2\*x\)\^\(1/4\)\) | Problem gg8 \(ix\) fourth roots of quadrat… | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_9\-i\-inst | PMT | x/\(\(2\*2\*x\^3\+3\)\*\(2\*x\^3\+3\)\^\(1/3\)\) | Problem gg9 \(i\) instantiated a=2,b=3 | FriCAS ELEM; book elem |
| generalised\_gunther\-prob\_gg\_9\-ii\-inst | PMT | x\^3/\(\(2\*2\*x\^3\+3\)\*\(2\*x\^3\+3\)\^\(1/3\)\) | Problem gg9 \(ii\) instantiated a=2,b=3 | FriCAS ELEM; book elem |
| generalised\_gunther\-prob\_gg\_9\-iii\_n1 | PMT | \(2\*a\*x\^3\+b\)\*x/\(\(a\^2\*m\*x\^6\+a\*b\*m\*x\^3\+b\^2\)\*\(a\*x\^3\+b\)\^\(1/3\)\) | Problem gg9 \(iii\) general form with n=1 | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_9\-iii\_n1\-ins… | PMT | \(2\*2\*x\^3\+3\)\*x/\(\(2\^2\*5\*x\^6\+2\*3\*5\*x\^3\+3\^2\)\*\(2\*x\^3\+3\)\^\(1/3\)\) | Problem gg9 \(iii\) general form with n=1 … | FriCAS ELEM; book elem |
| generalised\_gunther\-prob\_gg\_9\-iii\_n2 | PMT | \(3\*a\*x\^3\+2\*b\)\*x\^3/\(\(a\^2\*m\*x\^9\+a\*b\*m\*x\^6\+b\^2\)\*\(a\*x\^3\+b\)\^\(1/3\)\) | Problem gg9 \(iii\) general form with n=2 | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_9\-iii\_n2\-ins… | PMT | \(3\*2\*x\^3\+2\*3\)\*x\^3/\(\(2\^2\*5\*x\^9\+2\*3\*5\*x\^6\+3\^2\)\*\(2\*x\^3\+3\)\^\(1/3\)\) | Problem gg9 \(iii\) general form with n=2 … | FriCAS ELEM; book elem |
| generalised\_gunther\-prob\_gg\_9\-iii\_n3 | PMT | \(4\*a\*x\^3\+3\*b\)\*x\^5/\(\(a\^2\*m\*x\^12\+a\*b\*m\*x\^9\+b\^2\)\*\(a\*x\^3\+b\)\^\(1/3\)\) | Problem gg9 \(iii\) general form with n=3 | FriCAS ELEM; book elem |
| generalised\_gunther\-prob\_gg\_9\-iii\_n3\-ins… | PMT | \(4\*2\*x\^3\+3\*3\)\*x\^5/\(\(2\^2\*5\*x\^12\+2\*3\*5\*x\^9\+3\^2\)\*\(2\*x\^3\+3\)\^\(1/3\)\) | Problem gg9 \(iii\) general form with n=3 … | FriCAS ELEM; book elem |
| generalised\_gunther\-prob\_gg\_10\-n2 | PMT | 1/\(\(2\*b\*c/a\*x\^2\+c\)\*Sqrt\[a\+b\*x\^2\]\) | Problem gg10 Albert Rich form x\^\(n\-2\)/\(\(… | FriCAS ELEM; book elem |
| generalised\_gunther\-prob\_gg\_10\-n3\-inst | PMT | x/\(\(2\*3\*5/2\*x\^3\+5\)\*\(2\+3\*x\^3\)\^\(1/3\)\) | Problem gg10 Albert Rich form x\^\(n\-2\)/\(\(… | FriCAS ELEM; book elem |
| generalised\_gunther\-prob\_gg\_10\-n4 | PMT | x\^2/\(\(2\*b\*c/a\*x\^4\+c\)\*\(a\+b\*x\^4\)\^\(1/4\)\) | Problem gg10 Albert Rich form x\^\(n\-2\)/\(\(… | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_10\-n4\-inst | PMT | x\^2/\(\(2\*3\*5/2\*x\^4\+5\)\*\(2\+3\*x\^4\)\^\(1/4\)\) | Problem gg10 Albert Rich form x\^\(n\-2\)/\(\(… | FriCAS ELEM; book elem |
| generalised\_gunther\-prob\_gg\_11\-i\_plus | PMT | 1/\(\(x\+1\)\*\(x\^3\+x\)\^\(1/4\)\) | Problem gg11 \(i\) sign \+ | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_11\-i\_minus | PMT | 1/\(\(x\-1\)\*\(x\^3\+x\)\^\(1/4\)\) | Problem gg11 \(i\) sign \- | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_11\-ii\_minus | PMT | 1/\(\(x\-1\)\*\(x\^2\-2\*x\^3\)\^\(1/4\)\) | Problem gg11 \(ii\) signs \- \- | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_11\-iii\_plus | PMT | 1/\(\(a\*x\+b\)\*\(a\^2\*x\^3\+b\^2\*x\)\^\(1/4\)\) | Problem gg11 \(iii\) sign \+ | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_11\-iii\_plus\-… | PMT | 1/\(\(2\*x\+3\)\*\(2\^2\*x\^3\+3\^2\*x\)\^\(1/4\)\) | Problem gg11 \(iii\) sign \+ instantiated a… | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_11\-iii\_minus | PMT | 1/\(\(a\*x\-b\)\*\(a\^2\*x\^3\+b\^2\*x\)\^\(1/4\)\) | Problem gg11 \(iii\) sign \- | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_11\-iii\_minus… | PMT | 1/\(\(2\*x\-3\)\*\(2\^2\*x\^3\+3\^2\*x\)\^\(1/4\)\) | Problem gg11 \(iii\) sign \- instantiated a… | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_12\-i\_plus | PMT | 4\*x\^2/\(\(x\^4\+1\)\*\(x\^8\+1\)\^\(1/4\)\) | Problem gg12 \(i\) sign \+ \(variable u rena… | FriCAS ELEM; book elem |
| generalised\_gunther\-prob\_gg\_12\-i\_minus | PMT | 4\*x\^2/\(\(x\^4\-1\)\*\(x\^8\+1\)\^\(1/4\)\) | Problem gg12 \(i\) sign \- \(variable u rena… | FriCAS ELEM; book elem |
| generalised\_gunther\-prob\_gg\_12\-ii\_plus | PMT | 2/\(\(x\^2\+1\)\*\(1\+2\*x\^2\)\^\(1/4\)\) | Problem gg12 \(ii\) signs \+ \+ \(u renamed x… | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_12\-ii\_minus | PMT | 2/\(\(x\^2\-1\)\*\(1\-2\*x\^2\)\^\(1/4\)\) | Problem gg12 \(ii\) signs \- \- \(u renamed x… | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_12\-main | PMT | \(x\^3\+x\)\^\(1/4\)/\(\(x\-1\)\^2\*\(x\+1\)\) | Problem gg12 unnumbered target integral | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_12\-iii | PMT | 4\*x\^4\*\(x\^8\+1\)\^\(1/4\)/\(\(x\^4\-1\)\^2\*\(x\^4\+1\)\) | Problem gg12 \(iii\) \(u renamed x\) | FriCAS ELEM; book elem |
| generalised\_gunther\-prob\_gg\_13\-i\-inst | PMT | \(2\*5\*x\^4\+2\*3\*5\*x\^2\+3\*7\)\*Sqrt\[\(2\*x\^2\+3\)/\(5\*x\^2\+7\)\]/\(5\^2\*x\^8\+2\*5\*7\*x\^6\+\(7\^2\-2\^2\)\*x\^4\-2\*2\*3\*x… | Problem gg13 \(i\) rational radicand insta… | FriCAS ELEM; book elem |
| generalised\_gunther\-prob\_gg\_13\-ii | PMT | \(2\*a\*c\*x\^4\+\(a\*d\+3\*b\*c\)\*x\^2\+2\*b\*d\)\*\(\(a\*x\^2\+b\)/\(c\*x\^2\+d\)\)\^\(3/4\)/\(\(a\*x\^2\+b\)\*\(c\*x\^6\+d\*x\^4\-a\*x\^… | Problem gg13 \(ii\) rational radicand | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_13\-ii\-inst | PMT | \(2\*2\*5\*x\^4\+\(2\*7\+3\*3\*5\)\*x\^2\+2\*3\*7\)\*\(\(2\*x\^2\+3\)/\(5\*x\^2\+7\)\)\^\(3/4\)/\(\(2\*x\^2\+3\)\*\(5\*x\^6\+7\*x\^4\-2\*x\^… | Problem gg13 \(ii\) rational radicand inst… | FriCAS ELEM; book elem |
| generalised\_gunther\-prob\_gg\_13\-iii | PMT | \(2\*a\*c\*x\^3\+a\*d\*x\^2\+4\*b\*c\*x\+3\*b\*d\)/\(\(c\*x\+d\)\*\(c\*x\^4\+d\*x\^3\-a\*x\^2\-b\)\*\(\(a\*x\^2\+b\)/\(c\*x\+d\)\)\^\(1/3\)… | Problem gg13 \(iii\) rational radicand | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_13\-iii\-inst | PMT | \(2\*2\*5\*x\^3\+2\*7\*x\^2\+4\*3\*5\*x\+3\*3\*7\)/\(\(5\*x\+7\)\*\(5\*x\^4\+7\*x\^3\-2\*x\^2\-3\)\*\(\(2\*x\^2\+3\)/\(5\*x\+7\)\)\^\(1/3\)… | Problem gg13 \(iii\) rational radicand ins… | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_13\-iv | PMT | \(a\*c\*x\^6\+2\*b\*c\*x\^3\+b\*d\)\*\(\(a\*x\^3\+b\)/\(c\*x\^3\+d\)\)\^\(2/3\)/\(\(a\*x\^3\+b\)\*\(c\*x\^6\+\(d\-a\)\*x\^3\-b\)\) | Problem gg13 \(iv\) rational radicand | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_13\-iv\-inst | PMT | \(2\*5\*x\^6\+2\*3\*5\*x\^3\+3\*7\)\*\(\(2\*x\^3\+3\)/\(5\*x\^3\+7\)\)\^\(2/3\)/\(\(2\*x\^3\+3\)\*\(5\*x\^6\+\(7\-2\)\*x\^3\-3\)\) | Problem gg13 \(iv\) rational radicand inst… | FriCAS TIMEOUT; book elem |
| generalised\_gunther\-prob\_gg\_14\-i\_plus | PMT | 1/\(\(2\*x\+3\)\*\(1\-4\*x\^2\)\^\(1/3\)\) | Problem gg14 \(i\) sign \+ \(no statement\) | FriCAS TIMEOUT; book ? |
| generalised\_gunther\-prob\_gg\_14\-i\_minus | PMT | 1/\(\(2\*x\-3\)\*\(1\-4\*x\^2\)\^\(1/3\)\) | Problem gg14 \(i\) sign \- \(no statement\) | FriCAS TIMEOUT; book ? |
| generalised\_gunther\-prob\_gg\_14\-ii\_plus | PMT | 1/\(\(x\+1\)\*\(1\-9\*x\^2\)\^\(1/3\)\) | Problem gg14 \(ii\) sign \+ \(no statement\) | FriCAS TIMEOUT; book ? |
| generalised\_gunther\-prob\_gg\_14\-ii\_minus | PMT | 1/\(\(x\-1\)\*\(1\-9\*x\^2\)\^\(1/3\)\) | Problem gg14 \(ii\) sign \- \(no statement\) | FriCAS TIMEOUT; book ? |
| generalised\_gunther\-prob\_gg\_15\-ii\_plus | PMT | 1/\(\(2\*x\+1\)\*\(1\-8\*x\^3\)\^\(1/3\)\) | Problem gg15 \(ii\) 2x\+1 with 1\-8x\^3 \(no s… | FriCAS TIMEOUT; book ? |
| generalised\_gunther\-prob\_gg\_15\-ii\_minus | PMT | 1/\(\(2\*x\-1\)\*\(1\+8\*x\^3\)\^\(1/3\)\) | Problem gg15 \(ii\) 2x\-1 with 1\+8x\^3 \(no s… | FriCAS TIMEOUT; book ? |
| chebychev\-chebychev\_1\-iii | PMT | x\^\(4/3\)/\(a\*Sqrt\[x\]\+b\)\^\(2/3\) | Problem chebychev\_1 \(iii\) | FriCAS TIMEOUT; book elem |
| chebychev\-chebychev\_1\-iii\-inst | PMT | x\^\(4/3\)/\(2\*Sqrt\[x\]\+3\)\^\(2/3\) | Problem chebychev\_1 \(iii\) instantiated a… | FriCAS TIMEOUT; book elem |
| chebychev\-chebychev\_1\-iv | PMT | \(a\*x\^\(3/4\)\+b\)\^\(2/3\)/x\^3 | Problem chebychev\_1 \(iv\) | FriCAS TIMEOUT; book elem |
| chebychev\-chebychev\_1\-iv\-inst | PMT | \(2\*x\^\(3/4\)\+3\)\^\(2/3\)/x\^3 | Problem chebychev\_1 \(iv\) instantiated a=… | FriCAS TIMEOUT; book elem |
| chebychev\-chebychev\_3\-iii\-inst | PMT | x\^\(1/3\)\*\(2\*x\^2\+3\)\^\(1/3\) | Problem chebychev\_3 \(iii\) instantiated a… | FriCAS TIMEOUT; book elem |
| chebychev\-chebychev\_3\-v\-inst | PMT | 1/\(2\*x\^3\+3\)\^\(1/3\) | Problem chebychev\_3 \(v\) instantiated a=2… | FriCAS ELEM; book elem |
| chebychev\-chebychev\_3\-vi\-inst | PMT | x/\(2\*x\^3\+3\)\^\(1/3\) | Problem chebychev\_3 \(vi\) instantiated a=… | FriCAS NONELEM; book nonelem |
| chebychev\-chebychev\_3\-viii | PMT | x\^\(4/5\)\*\(a/x\^\(4/5\)\+b\)\^\(5/4\) | Problem chebychev\_3 \(viii\) \(expect by bo… | FriCAS TIMEOUT; book elem |
| chebychev\-chebychev\_3\-x | PMT | x\^\(3/2\)/Sqrt\[a/Sqrt\[x\]\+b\] | Problem chebychev\_3 \(x\) \(expect by book … | FriCAS ELEM; book elem |
| chebychev\-chebychev\_3c\-i\_n3 | PMT | 1/\(Sin\[x\]\)\^\(1/3\) | Problem chebychev\_3c \(i\) 1/sin\^\(1/n\) n=3… | FriCAS NONELEM; book nonelem |
| chebychev\-chebychev\_3c\-iv\_n3 | PMT | \(Sin\[x\]\)\^\(1/3\) | Problem chebychev\_3c \(iv\) 1/csc\^\(1/n\)=si… | FriCAS NONELEM; book nonelem |
| chebychev\-chebychev\_4\-iv\_n3 | PMT | \(1/\(\(Exp\[x\]\+Exp\[\-x\]\)/2\)\)\^\(1/3\) | Problem chebychev\_4 \(iv\) sech\^\(1/n\) via … | FriCAS NONELEM; book nonelem |
| chebychev\-chebychev\_5\-v | PMT | \(a\*Sqrt\[x\]\+b\)\^\(2/3\)/\(c\*Sqrt\[x\]\+d\)\^\(2/3\) | Problem chebychev\_5 \(v\) type IV conjectu… | FriCAS TIMEOUT; book elem |
| chebychev\-chebychev\_5\-v\-inst | PMT | \(2\*Sqrt\[x\]\+3\)\^\(2/3\)/\(5\*Sqrt\[x\]\+7\)\^\(2/3\) | Problem chebychev\_5 \(v\) type IV conjectu… | FriCAS TIMEOUT; book elem |
| chebychev\-chebychev\_5\-viii | PMT | \(a\*x\^\(1/3\)\+b\)\^\(5/3\)/\(c\*x\^\(1/3\)\+d\)\^\(2/3\) | Problem chebychev\_5 \(viii\) type V conjec… | FriCAS TIMEOUT; book elem |
| chebychev\-chebychev\_5\-viii\-inst | PMT | \(2\*x\^\(1/3\)\+3\)\^\(5/3\)/\(5\*x\^\(1/3\)\+7\)\^\(2/3\) | Problem chebychev\_5 \(viii\) type V conjec… | FriCAS TIMEOUT; book elem |
| inverse\_function\-inverse\_function\_2\-v | PMT | \(2\^\(1/3\)\*\(Sqrt\[x\^2\+4\]\+x\)\^\(2/3\)\-2\)\^\(1/3\)/\(Sqrt\[x\^2\+4\]\+x\)\^\(1/9\) | Problem inverse\_function\_2 \(v\) | FriCAS TIMEOUT; book elem |
| realis\-prob\_realis\_4\-ii | PMT | \(a\-2\*x\^2\)/\(\(2\*x\^4\-2\*a\*x\^2\+a\)\*Sqrt\[\(a\-x\^2\)\*\(x\^4\-a\*x\^2\+a\)\]\) | Problem prob\_realis\_4 \(ii\) \(a=2 gives \(i… | FriCAS ELEM; book elem |
| nested\-nested\_3\-i\-inst | PMT | \(2/x\^3\+3\*x\*Sqrt\[x\]\)\^\(1/3\) | Problem nested\_3 \(i\), n=\-3; inst a=2, b=… | FriCAS TIMEOUT; book elem |
| nested\-nested\_3\-iii | PMT | \(a\*x\+b\*x\*Sqrt\[x\]\)\^\(1/3\) | Problem nested\_3 \(iii\), n=1 | FriCAS TIMEOUT; book elem |
| nested\-nested\_3\-iii\-inst | PMT | \(2\*x\+3\*x\*Sqrt\[x\]\)\^\(1/3\) | Problem nested\_3 \(iii\), n=1; inst a=2, b… | FriCAS TIMEOUT; book elem |
| nested\-nested\_3\-v\-inst | PMT | \(2\*x\^3\+3\*x\*Sqrt\[x\]\)\^\(1/3\) | Problem nested\_3 \(v\), n=3; inst a=2, b=3 | FriCAS ELEM; book elem |
| nested\-nested\_4\-ex1 | PMT | \(a/x\^4\+b\*x\^\(4/3\)\)\^\(1/4\) | Problem nested\_4, solution example 1 \(n=… | FriCAS TIMEOUT; book elem |
| nested\-nested\_4\-ex1\-inst | PMT | \(2/x\^4\+3\*x\^\(4/3\)\)\^\(1/4\) | Problem nested\_4, solution example 1 \(n=… | FriCAS TIMEOUT; book elem |
| nested\-nested\_6\-i | PMT | Sqrt\[d\*x\^2\+c\*x\^2/Sqrt\[a/x\+b\]\]/Sqrt\[a/x\+b\] | Problem nested\_6 \(i\) | FriCAS ELEM; book elem |
| nested\-nested\_6\-ii | PMT | Sqrt\[a/x\+b\]/Sqrt\[d\*x\^2\+c\*x\^2/Sqrt\[a/x\+b\]\] | Problem nested\_6 \(ii\) | FriCAS ELEM; book elem |
| nested\-nested\_6\-iii | PMT | Sqrt\[a/x\+b\]\*Sqrt\[d\*x\^2\+c\*x\^2/Sqrt\[a/x\+b\]\] | Problem nested\_6 \(iii\) | FriCAS ELEM; book elem |
| nested\-nested\_6\-iv | PMT | 1/\(Sqrt\[a/x\+b\]\*Sqrt\[d\*x\^2\+c\*x\^2/Sqrt\[a/x\+b\]\]\) | Problem nested\_6 \(iv\) | FriCAS ELEM; book elem |
| nested\-nested\_9\-i\-inst | PMT | 1/\(\-2\*x\+Sqrt\[2\^2\*x\^2\+3\*x\+5\]\)\^\(1/4\) | Problem nested\_9 \(i\); inst a=2, b=3, c=5 | FriCAS ELEM; book elem |
| nested\-nested\_9\-ii\-inst | PMT | 1/\(\-2\*x\+Sqrt\[2\^2\*x\^2\+3\*x\+5\]\)\^\(3/4\) | Problem nested\_9 \(ii\); inst a=2, b=3, c=… | FriCAS ELEM; book elem |
| nested\-nested\_9\-iii\-inst | PMT | \(\-2\*x\+Sqrt\[2\^2\*x\^2\+3\*x\+5\]\)\^\(1/4\)/\(3\*x\+5\) | Problem nested\_9 \(iii\); inst a=2, b=3, c… | FriCAS ELEM; book elem |
| nested\-nested\_9\-iv\-inst | PMT | 1/\(\-2\*x\+Sqrt\[2\^2\*x\^2\+3\*x\+5\]\)\^\(5/6\) | Problem nested\_9 \(iv\); inst a=2, b=3, c=… | FriCAS ELEM; book elem |
| nested\-nested\_9\-v\-inst | PMT | x/\(\-2\*x\+Sqrt\[2\^2\*x\^2\+3\*x\+5\]\)\^\(1/8\) | Problem nested\_9 \(v\); inst a=2, b=3, c=5 | FriCAS ELEM; book elem |
| nested\-nested\_12\-ii\-inst | PMT | \(2\*x\-Sqrt\[2\^2\*x\^2\+3\]\)\^\(1/3\)/\(x\*Sqrt\[2\^2\*x\^2\+3\]\) | Problem nested\_12 \(ii\); inst a=2, b=3 | FriCAS ELEM; book elem |
| nested\-nested\_12\-iii\-inst | PMT | \(2\*x\-Sqrt\[2\^2\*x\^2\+3\]\)\^\(1/3\)/\(\(2\^2\*x\^2\+3\)\^2\*Sqrt\[2\^2\*x\^2\+3\]\) | Problem nested\_12 \(iii\); inst a=2, b=3 | FriCAS ELEM; book elem |
| nested\-nested\_14\-iv\-inst | PMT | Sqrt\[Sqrt\[2\*x\^2\+3\^2\]\+3\]/\(\(2\*x\^2\-3\)\*Sqrt\[2\*x\^2\+3\^2\]\) | Problem nested\_14 \(iv\); inst a=2, b=3 | FriCAS TIMEOUT; book elem |
| nested\-nested\_14\-v\-inst | PMT | \(2\*x\^2\+3\)\^2\*Sqrt\[2\*x\^2\+3\^2\]/\(\(2\*x\^2\-3\)\^2\*Sqrt\[Sqrt\[2\*x\^2\+3\^2\]\+3\]\) | Problem nested\_14 \(v\); inst a=2, b=3 | FriCAS TIMEOUT; book elem |
| nested\-nested\_17\-i | PMT | \(\(x\^2\+1\)\^\(3/2\)\+x\^3\)\^\(1/3\)/\(x\^2\*\(x\^2\+1\)\) | Problem nested\_17 \(i\) | FriCAS TIMEOUT; book elem |
| nested\-nested\_17\-iii | PMT | \(2\*x\^2\+1\)\*\(x\^3\*\(x\^2\+1\)\^\(3/2\)\+1\)\^\(1/3\)/\(x\*\(x\^2\+1\)\) | Problem nested\_17 \(iii\) | FriCAS ERROR; book elem |
| nested\-nested\_19\-xi | PMT | \(2\+1/x\^2\-Sqrt\[8\*x\^2\+1\]/x\^2\)\^\(1/4\) | Problem nested\_19 \(xi\) | FriCAS TIMEOUT; book elem |
| multrad\-multiple\_radicals\_0\-v | PMT | x\^2\*\(x\^4\-2\*b\*x\^2\+a\*b\)\*Sqrt\[x\^2\-a\]/\(\(d\*x\^8\-2\*a\*d\*x\^6\+\(a\^2\*d\+c\)\*x\^4\-2\*b\*c\*x\^2\+b\^2\*c\)\*Sqrt\[x\^… | Problem multiple\_radicals\_0 \(v\) | FriCAS TIMEOUT; book elem |
| multrad\-multiple\_radicals\_0\-vi | PMT | \(x\^4\+2\*a\*x\^2\+a\^2\)\*\(x\^6\+3\*a\*x\^4\-a\*\(2\*a\+3\*b\)\*x\^2\+a\^2\*b\)/\(\(c\*x\^12\-2\*c\*x\^10\+\(a\^2\*c\+4\*a\*b\*c\+b\^2… | Problem multiple\_radicals\_0 \(vi\) | FriCAS TIMEOUT; book elem |
| multrad\-multiple\_radicals\_0\-vi\-inst | PMT | \(x\^4\+2\*2\*x\^2\+2\^2\)\*\(x\^6\+3\*2\*x\^4\-2\*\(2\*2\+3\*3\)\*x\^2\+2\^2\*3\)/\(\(5\*x\^12\-2\*5\*x\^10\+\(2\^2\*5\+4\*2\*3\*5\+3\^2… | Problem multiple\_radicals\_0 \(vi\); inst a… | FriCAS TIMEOUT; book elem |
| multrad\-multiple\_radicals\_2\-iii | PMT | \(x\^4\+2\*a\*x\^2\-a\*b\)\*Sqrt\[x\^2\+a\]/\(\(x\^4\-\(b\+1\)\*x\^2\-a\)\*\(x\^4\+\(1\-b\)\*x\^2\+a\)\*Sqrt\[x\^2\-b\]\) | Problem multiple\_radicals\_2 \(iii\) | FriCAS ELEM; book ? |
| multrad\-multiple\_radicals\_3\-i | PMT | \(x\^2\-2\*a\*x\+a\*b\+a\*c\-b\*c\)\*Sqrt\[x\-b\]\*Sqrt\[x\-c\]/\(\(x\^2\-\(b\+c\-1\)\*x\+b\*c\-a\)\*\(x\-a\)\^\(3/2\)\) | Problem multiple\_radicals\_3 \(i\) | FriCAS TIMEOUT; book ? |
| multrad\-multiple\_radicals\_3\-iii | PMT | \(x\^2\-2\*b\*x\+a\*b\-a\*c\+b\*c\)\*Sqrt\[x\-b\]/\(\(x\^2\-\(a\+c\-1\)\*x\-b\+a\*c\)\*\(x\^2\-\(a\+c\+1\)\*x\+b\+a\*c\)\*Sqrt\[x\-a\]\*S… | Problem multiple\_radicals\_3 \(iii\) | FriCAS ELEM; book ? |
| multrad\-multiple\_radicals\_7\-iv | PMT | ArcTan\[d\*Sqrt\[a\*x\+b\]\-d\*Sqrt\[a\*x\+c\]\] | Problem multiple\_radicals\_7 \(iv\) | FriCAS ELEM; book elem |
| euler\-prob\_euler\_outrageous\_1\-ii | PMT | \(a\^2\*x\^2\-2\*a\*b\*x\+b\^2\)\*\(a\^2\*x\^2\+2\*a\*b\*x\+b\^2\)/\(\(a\^2\*x\^2\+b\^2\)\*\(a\^4\*x\^4\+6\*a\^2\*b\^2\*x\^2\+b\^4\)\^\(3/… | Problem prob\_euler\_outrageous\_1 \(ii\) | FriCAS TIMEOUT; book elem |
| euler\-prob\_euler\_outrageous\_1\-ii\-inst | PMT | \(2\^2\*x\^2\-2\*2\*3\*x\+3\^2\)\*\(2\^2\*x\^2\+2\*2\*3\*x\+3\^2\)/\(\(2\^2\*x\^2\+3\^2\)\*\(2\^4\*x\^4\+6\*2\^2\*3\^2\*x\^2\+3\^4\)\^\(3/… | Problem prob\_euler\_outrageous\_1 \(ii\); in… | FriCAS TIMEOUT; book elem |
| euler\-prob\_euler\_outrageous\_2\-m | PMT | 1/\(\(3\-x\^2\)\*\(1\-3\*x\^2\)\^\(1/3\)\) | Problem prob\_euler\_outrageous\_2 \(\- sign\)… | FriCAS TIMEOUT; book elem |
| euler\-prob\_euler\_outrageous\_3\-ii\-inst | Trager | \(2\*x\^2\-x\-2\*3\)/\(\(x\^4\-2\*x\^2\+3\)\*\(2\*x\^2\-3\)\^\(1/4\)\) | Problem prob\_euler\_outrageous\_3 \(ii\); in… | FriCAS TIMEOUT; book elem |
| euler\-prob\_euler\_outrageous\_4\-iii | PMT | Sqrt\[a\*x\^4\+c\*x\^2\+b\]/\(b\-a\*x\^4\) | Problem prob\_euler\_outrageous\_4 \(iii\) | FriCAS ELEM; book elem |
| euler\-prob\_euler\_outrageous\_4\-iv | PMT | x\^2/\(\(b\-a\*x\^4\)\*Sqrt\[a\*x\^4\+b\]\) | Problem prob\_euler\_outrageous\_4 \(iv\) | FriCAS ELEM; book elem |
| euler\-prob\_euler\_outrageous\_4\-iv\-inst | PMT | x\^2/\(\(3\-2\*x\^4\)\*Sqrt\[2\*x\^4\+3\]\) | Problem prob\_euler\_outrageous\_4 \(iv\); in… | FriCAS ELEM; book elem |
| euler\-prob\_euler\_outrageous\_5 | PMT | \(a\*x\^2\-b\)/\(\(a\*x\^2\+b\)\*Sqrt\[a\^2\*x\^4\+b\^2\]\) | Problem prob\_euler\_outrageous\_5 | FriCAS ELEM; book elem |
| timofeev\-timofeev\_2\-i | PMT | 1/\(\(x\^3\-1\)\*\(x\^3\+2\)\^\(1/3\)\) | Problem timofeev\_2 \(i\) | FriCAS ELEM; book elem |
| timofeev\-timofeev\_2\-v | PMT | \(x\^5\-2\)\^2/\(x\^5\+3\)\^\(16/5\) | Problem timofeev\_2 \(v\) | FriCAS ELEM; book elem |
| timofeev\-timofeev\_2\-vi | PMT | \(x\^5\-2\)\^2/\(\(x\^5\+3\)\^3\*\(x\^5\+3\)\^\(1/5\)\) | Problem timofeev\_2 \(vi\) | FriCAS ELEM; book elem |
| timofeev\-timofeev\_2\-vii | PMT | 1/\(\(x\^3\+3\*x\^2\+3\*x\)\*\(x\^3\+3\*x\^2\+3\*x\+3\)\^\(1/3\)\) | Problem timofeev\_2 \(vii\) | FriCAS ELEM; book elem |
| timofeev\-timofeev\_3\-v | PMT | \(\(x\^2\+2\)/x\^2\)\^\(7/9\)/\(x\^2\+2\)\^\(3/2\) | Problem timofeev\_3 \(v\) | FriCAS ELEM; book ? |
| potpourri\-p6\-v | PMT | 1/\(\(x\-1\)\*\(x\^3\+x\)\^\(1/4\)\) | Potpourri problem 6 \(v\) | FriCAS TIMEOUT; book ? |
| potpourri\-p6\-xi | PMT | \(x\^4\+4\)\*\(x\^4\-1\)\^\(3/4\)/\(x\^8\*\(x\^4\-4\)\) | Potpourri problem 6 \(xi\) | FriCAS ELEM; book ? |
| potpourri\-p6\-xiii | PMT | \(x\-2\)/\(\(x\^2\+2\)\*\(2\*x\^2\+x\-1\)\^\(1/3\)\) | Potpourri problem 6 \(xiii\) | FriCAS TIMEOUT; book ? |
| potpourri\-p6\-xiv | PMT | \(x\^3\+2\)\*\(2\*x\^3\+1\)\^\(2/3\)/\(x\^6\*\(x\^3\+1\)\) | Potpourri problem 6 \(xiv\) | FriCAS TIMEOUT; book ? |
| potpourri\-p6\-xv | PMT | \(x\^6\+x\^3\+1\)/\(\(1\-x\^6\)\*\(x\^5\+x\^3\)\^\(1/4\)\) | Potpourri problem 6 \(xv\) | FriCAS ERROR; book ? |
| potpourri\-p6\-xvi | PMT | x\^3/\(\(x\^6\-1\)\*\(x\^4\+x\^2\)\^\(1/3\)\) | Potpourri problem 6 \(xvi\) | FriCAS TIMEOUT; book ? |
| danger\-d1\-i\-p | PMT | 1/\(\(3\*x\^4\-4\)\*\(x\^2\+1\)\^\(1/4\)\) | Danger zone problem 1 \(i\), x\^2\+1; book t… | FriCAS TIMEOUT; book elem |
| danger\-d1\-i\-m | PMT | 1/\(\(3\*x\^4\-4\)\*\(x\^2\-1\)\^\(1/4\)\) | Danger zone problem 1 \(i\), x\^2\-1; book t… | FriCAS TIMEOUT; book elem |
| danger\-d1\-ii\-p | PMT | x\^2/\(\(3\*x\^4\-4\)\*\(x\^2\+1\)\^\(1/4\)\) | Danger zone problem 1 \(ii\), x\^2\+1; book … | FriCAS TIMEOUT; book elem |
| danger\-d1\-ii\-m | PMT | x\^2/\(\(3\*x\^4\-4\)\*\(x\^2\-1\)\^\(1/4\)\) | Danger zone problem 1 \(ii\), x\^2\-1; book … | FriCAS TIMEOUT; book elem |
| danger\-d1\-iii\-p | PMT | 1/\(\(x\^4\+4\)\*\(x\^2\+1\)\^\(1/4\)\) | Danger zone problem 1 \(iii\), x\^2\+1; book… | FriCAS TIMEOUT; book elem |
| danger\-d1\-iii\-m | PMT | 1/\(\(x\^4\+4\)\*\(x\^2\-1\)\^\(1/4\)\) | Danger zone problem 1 \(iii\), x\^2\-1; book… | FriCAS TIMEOUT; book elem |
| danger\-d1\-iv\-p | PMT | x\^2/\(\(x\^4\+4\)\*\(x\^2\+1\)\^\(1/4\)\) | Danger zone problem 1 \(iv\), x\^2\+1; book … | FriCAS TIMEOUT; book elem |
| danger\-d1\-iv\-m | PMT | x\^2/\(\(x\^4\+4\)\*\(x\^2\-1\)\^\(1/4\)\) | Danger zone problem 1 \(iv\), x\^2\-1; book … | FriCAS TIMEOUT; book elem |
| danger\-d1\-v\-p | Trager | 1/\(\(x\+1\)\*\(1\-5\*x\^2\)\^\(1/5\)\) | Danger zone problem 1 \(v\), x\+1 | FriCAS ERROR; book elem |
| danger\-d1\-v\-p | PMT | 1/\(\(x\+1\)\*\(1\-5\*x\^2\)\^\(1/5\)\) | Danger zone problem 1 \(v\), x\+1 | FriCAS ERROR; book elem |
| danger\-d1\-v\-m | Trager | 1/\(\(x\-1\)\*\(1\-5\*x\^2\)\^\(1/5\)\) | Danger zone problem 1 \(v\), x\-1 | FriCAS ERROR; book elem |
| danger\-d1\-v\-m | PMT | 1/\(\(x\-1\)\*\(1\-5\*x\^2\)\^\(1/5\)\) | Danger zone problem 1 \(v\), x\-1 | FriCAS ERROR; book elem |
| multmeth\-m1 | PMT | \(a\^3\*x\^3\-b\^3\)\^\(1/3\)/\(a\*x\+b\) | Multiple methods problem 1 \[4\] | FriCAS ERROR; book ? |
| multmeth\-m1\-inst | PMT | \(2\^3\*x\^3\-3\^3\)\^\(1/3\)/\(2\*x\+3\) | Multiple methods problem 1 \[4\]; inst a=2… | FriCAS TIMEOUT; book ? |
| multmeth\-m2\-inst | PMT | \(2\*x\+3\)/\(\(2\-x\^2\)\*\(x\^2\-1\)\^\(1/4\)\) | Multiple methods problem 2; inst a=2, b=… | FriCAS TIMEOUT; book ? |
| intro\-multiple\_approaches\-1 | PMT | 1/\(\(x\-1\)\*Sqrt\[x\^2\+b\*x\+c\]\) | Introduction, Multiple Approaches to Int… | FriCAS ELEM; book elem |

### Elementary per book or FriCAS, but declined (606)

| id | method | integrand | book ref | detail |
|-----------|------|------------------|---------|--------------|
| warmup\-prob\_warmup\_2b\-plus\-n2 | PMT | \(\(2\*2\+1\)\*x\-1\)/\(\(x\^2\-1\)\*Sqrt\[\(x\+1\)\*\(x\^2\-1\)\^2\-1\]\) | Problem 0.2b Realis upper sign n=2 | \{"needs torsion realisation \(milestone iii\)", \{\{\-1 \+ Sqrt\[x\], \{\-2\*I, 2\*I\}\}, \{1 \+ Sqrt\[x\], \{2\*I, \-2\*I\}\}, \{1 \+ x… |
| warmup\-prob\_warmup\_2b\-minus\-n2 | PMT | \(\(2\*2\+1\)\*x\+1\)/\(\(x\^2\-1\)\*Sqrt\[\(x\-1\)\*\(x\^2\-1\)\^2\-1\]\) | Problem 0.2b Realis lower sign n=2 | \{"needs torsion realisation \(milestone iii\)", \{\{\-1 \+ x, \{\-3\*I, 3\*I\}\}, \{1 \+ x, \{\-2\*I, 2\*I\}\}\}\} |
| warmup\-prob\_warmup\_2b\-plus\-n3 | PMT | \(\(2\*3\+1\)\*x\-1\)/\(\(x\^2\-1\)\*Sqrt\[\(x\+1\)\*\(x\^2\-1\)\^3\-1\]\) | Problem 0.2b Realis upper sign n=3 | \{"needs torsion realisation \(milestone iii\)", \{\{\-1 \+ x, \{\-3\*I, 3\*I\}\}, \{1 \+ x, \{\-4\*I, 4\*I\}\}\}\} |
| warmup\-prob\_warmup\_2b\-minus\-n3 | PMT | \(\(2\*3\+1\)\*x\+1\)/\(\(x\^2\-1\)\*Sqrt\[\(x\-1\)\*\(x\^2\-1\)\^3\-1\]\) | Problem 0.2b Realis lower sign n=3 | \{"needs torsion realisation \(milestone iii\)", \{\{\-1 \+ Sqrt\[x\], \{\-4\*I, 4\*I\}\}, \{1 \+ Sqrt\[x\], \{4\*I, \-4\*I\}\}, \{1 \+ x… |
| warmup\-prob\_warmup\_2c\-i | Trager | 1/\(1\-x\^3\) | Problem 0.2c \(i\) rational | \{"failed", "no radical in the integrand \(a rational function\)"\} |
| warmup\-prob\_warmup\_2c\-ii | Trager | 1/\(x\^4\+1\) | Problem 0.2c \(ii\) rational | \{"failed", "no radical in the integrand \(a rational function\)"\} |
| warmup\-prob\_warmup\_2c\-iii | Trager | x\^4/\(x\^6\-b\^6\) | Problem 0.2c \(iii\) rational | \{"failed", "integrand is not algebraic over Q\(x\)"\} |
| warmup\-prob\_warmup\_2c\-iii\-inst | Trager | x\^4/\(x\^6\-3\^6\) | Problem 0.2c \(iii\) rational instantiated… | \{"failed", "no radical in the integrand \(a rational function\)"\} |
| warmup\-prob\_warmup\_2c\-iv | Trager | \(x\^4\+a\)/\(x\^4\+b\)\^2 | Problem 0.2c \(iv\) rational | \{"failed", "integrand is not algebraic over Q\(x\)"\} |
| warmup\-prob\_warmup\_2c\-iv\-inst | Trager | \(x\^4\+2\)/\(x\^4\+3\)\^2 | Problem 0.2c \(iv\) rational instantiated … | \{"failed", "no radical in the integrand \(a rational function\)"\} |
| warmup\-prob\_warmup\_2cc\-i | Trager | x/\(x\^8\+1\) | Problem 1.2cc \(i\) rational | \{"failed", "no radical in the integrand \(a rational function\)"\} |
| warmup\-prob\_warmup\_2cc\-i | PMT | x/\(x\^8\+1\) | Problem 1.2cc \(i\) rational | \{"failed", "residue class lost: no places over the prime for a root of the residue polynomial", 1 \+ x\^8, \(\-1/8… |
| warmup\-prob\_warmup\_2cc\-ii | Trager | x\^2/\(x\^4\+x\^2\+1\) | Problem 1.2cc \(ii\) rational | \{"failed", "no radical in the integrand \(a rational function\)"\} |
| warmup\-prob\_warmup\_2cc\-iii | Trager | \(x\^2\+2\*x\)/\(x\^4\+2\*x\^2\+4\*x\+2\) | Problem 1.2cc \(iii\) rational | \{"failed", "no radical in the integrand \(a rational function\)"\} |
| warmup\-prob\_warmup\_2cc\-iv | Trager | x\^2/\(x\^6\-3\*x\^4\+4\*x\^2\-2\) | Problem 1.2cc \(iv\) rational | \{"failed", "no radical in the integrand \(a rational function\)"\} |
| warmup\-prob\_warmup\_2d | Trager | \(x\^6\+\(3\-2\*Sqrt\[2\]\)\*x\^4\+\(12\*Sqrt\[2\]\-17\)\*x\^2\+70\*Sqrt\[2\]\-99\)/\(x\^8\+\(8\-8\*Sqrt\[2\]\)\*x\^6\+\(78\-56\*Sq… | Problem 1.2d rational over Q\(Sqrt\[2\]\) | \{"failed", "algebraic constant in the integrand \(coefficients must be rational\)"\} |
| warmup\-prob\_warmup\_2e | Trager | \(189\*x\^5\+945\*x\^4\+1131\*x\^3\-147\*x\^2\-110\*x\-28\)/\(9\*x\^6\+54\*x\^5\+74\*x\^4\-14\*x\^3\+25\*x\^2\-8\*x\+4\) | Problem 2.2e \(and 2.2f\) Lazard\-Rioboo\-Tr… | \{"failed", "no radical in the integrand \(a rational function\)"\} |
| warmup\-prob\_warmup\_2e | PMT | \(189\*x\^5\+945\*x\^4\+1131\*x\^3\-147\*x\^2\-110\*x\-28\)/\(9\*x\^6\+54\*x\^5\+74\*x\^4\-14\*x\^3\+25\*x\^2\-8\*x\+4\) | Problem 2.2e \(and 2.2f\) Lazard\-Rioboo\-Tr… | \{"failed", "residue class lost: no places over the prime for a root of the residue polynomial", 4 \- 8 x \+ 25 x… |
| warmup\-prob\_warmup\_2g | Trager | \(\(1\+Sqrt\[3\]\)\*x\^2\-4\*Sqrt\[3\]\+4\)/\(\(9\+5\*Sqrt\[3\]\)\*x\^4\+24\*x\^2\-80\*Sqrt\[3\]\+144\) | Problem 1.2g rational over Q\(Sqrt\[3\]\) | \{"failed", "algebraic constant in the integrand \(coefficients must be rational\)"\} |
| warmup\-prob\_warmup\_3\-i | Trager | \(x\^2\-a\)/\(Sqrt\[x\+1\]\+b\) | Problem 0.3 \(i\) Q\(x,\(ax\+b\)\^\(m/n\)\) | \{"failed", "integrand is not a rational function of x and the radical over Q"\} |
| warmup\-prob\_warmup\_3\-ii | Trager | \(\(a\*x\+b\)\^\(1/3\)\+c\)/\(\(a\*x\+b\)\^\(1/3\)\+d\) | Problem 0.3 \(ii\) Q\(x,\(ax\+b\)\^\(m/n\)\) | \{"failed", "radicand is not a rational function over Q"\} |
| warmup\-prob\_warmup\_3\-iii | Trager | \(Sqrt\[a\*x\+b\]\+a\)/\(\(a\*x\+b\)\^\(1/4\)\+b\) | Problem 0.3 \(iii\) Q\(x,\(ax\+b\)\^\(m/n\)\) | \{"failed", "radicand is not a rational function over Q"\} |
| warmup\-prob\_warmup\_4\-i | Trager | 1/Sqrt\[a\*x\^2\+b\*x\+c\] | Problem 0.4 \(i\) Euler substitution | \{"failed", "radicand is not a rational function over Q"\} |
| warmup\-prob\_warmup\_4\-i | PMT | 1/Sqrt\[a\*x\^2\+b\*x\+c\] | Problem 0.4 \(i\) Euler substitution | \{"failed", "no solution within bounds", \{0\}\} |
| warmup\-prob\_warmup\_4\-ii | Trager | Sqrt\[a\*x\^2\+b\*x\+c\] | Problem 0.4 \(ii\) Euler substitution | \{"failed", "radicand is not a rational function over Q"\} |
| warmup\-prob\_warmup\_4\-ii | PMT | Sqrt\[a\*x\^2\+b\*x\+c\] | Problem 0.4 \(ii\) Euler substitution | \{"failed", "no solution within bounds", \{2\}\} |
| warmup\-prob\_warmup\_4\-iii | Trager | Sqrt\[a\*x\^2\+b\*x\+c\]/x\^2 | Problem 0.4 \(iii\) Euler substitution | \{"failed", "radicand is not a rational function over Q"\} |
| warmup\-prob\_warmup\_4\-iv | Trager | x\^2\*Sqrt\[a\*x\^2\+b\*x\+c\] | Problem 0.4 \(iv\) Euler substitution | \{"failed", "radicand is not a rational function over Q"\} |
| warmup\-prob\_warmup\_4\-v | Trager | \(a\*x\^2\+b\*x\+c\)\^\(3/2\)/x\^2 | Problem 0.4 \(v\) Euler substitution | \{"failed", "radicand is not a rational function over Q"\} |
| warmup\-prob\_warmup\_4\-vi | Trager | x\*\(a\*x\^2\+b\*x\+c\)\^\(3/2\) | Problem 0.4 \(vi\) Euler substitution | \{"failed", "radicand is not a rational function over Q"\} |
| warmup\-prob\_warmup\_4\-vii | Trager | \(a\*x\^2\+b\*x\+c\)\^\(5/2\)/x | Problem 0.4 \(vii\) Euler substitution | \{"failed", "radicand is not a rational function over Q"\} |
| warmup\-prob\_warmup\_4\-viii | Trager | x\*\(a\*x\^2\+b\*x\+c\)\^\(5/2\) | Problem 0.4 \(viii\) Euler substitution | \{"failed", "radicand is not a rational function over Q"\} |
| warmup\-prob\_warmup\_4\-ix | Trager | 1/\(\(x\+b\)\^2\*\(a\*x\^2\+b\*x\+c\)\^\(3/2\)\) | Problem 0.4 \(ix\) Euler substitution | \{"failed", "radicand is not a rational function over Q"\} |
| warmup\-prob\_warmup\_4\-x | Trager | 1/\(\(x\+b\)\^2\*\(a\*x\^2\+b\*x\+c\)\^\(5/2\)\) | Problem 0.4 \(x\) Euler substitution | \{"failed", "radicand is not a rational function over Q"\} |
| warmup\-prob\_warmup\_4\-xi | Trager | x/\(\(a\*x\^2\+b\*x\+c\)\^\(3/2\)\-\(a\*x\^2\+b\*x\+c\)\^\(5/2\)\) | Problem 0.4 \(xi\) Euler substitution | \{"failed", "radicand is not a rational function over Q"\} |
| warmup\-prob\_warmup\_5 | Trager | \(\(a\*x\+b\)\^\(2/3\)\+1\)/\(\(a\*x\+b\)\^\(1/2\)\+\(a\*x\+b\)\^\(1/4\)\) | Problem 0.5 multiple radicals of ax\+b | \{"failed", "radicand is not a rational function over Q"\} |
| warmup\-prob\_warmup\_7 | Trager | \(\(x\-5\)\*\(3\*x\-1\)\^\(5/2\)\)/\(\(x\+2\)\*\(8\-7\*x\)\^\(3/2\)\) | Problem 0.7 two linear square roots | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| warmup\-prob\_warmup\_7 | PMT | \(\(x\-5\)\*\(3\*x\-1\)\^\(5/2\)\)/\(\(x\+2\)\*\(8\-7\*x\)\^\(3/2\)\) | Problem 0.7 two linear square roots | \{"failed", "no solution within bounds", \{5\}\} |
| warmup\-prob\_warmup\_10\-i | Trager | 1/\(1\-Sin\[x\]\) | Problem 0.10 \(i\) tangent half\-angle | \{"failed", "transcendental function in the integrand"\} |
| warmup\-prob\_warmup\_10\-ii | Trager | \(Sin\[x\]\-1\)/\(Sin\[x\]\+Cos\[x\]\) | Problem 0.10 \(ii\) tangent half\-angle | \{"failed", "transcendental function in the integrand"\} |
| warmup\-prob\_warmup\_10\-iii | Trager | \(1/Cos\[x\]\-Cos\[x\]/Sin\[x\]\)/\(Cos\[x\]\-Sin\[x\]\+1\) | Problem 0.10 \(iii\) tangent half\-angle, s… | \{"failed", "transcendental function in the integrand"\} |
| warmup\-prob\_warmup\_13\-iii | Trager | 1/\(x\*Sqrt\[a\*x\^2\+b\*x\+c\]\) | Problem 1.13 \(iii\) | \{"failed", "radicand is not a rational function over Q"\} |
| warmup\-prob\_warmup\_14\-ii | Trager | \(\(Exp\[x\]\-Exp\[\-x\]\)/2\)/\(1\+\(Exp\[x\]\-Exp\[\-x\]\)/2\) | Problem 1.14 \(ii\) sinh\(x\)/\(1\+sinh\(x\)\), s… | \{"failed", "transcendental function in the integrand"\} |
| warmup\-prob\_warmup\_14\-iii | Trager | 1/\(\(\(Exp\[x\]\+Exp\[\-x\]\)/2\)\^3\-\(\(Exp\[x\]\+Exp\[\-x\]\)/2\)\) | Problem 1.14 \(iii\) 1/\(cosh\(x\)\^3\-cosh\(x\)\)… | \{"failed", "transcendental function in the integrand"\} |
| warmup\-prob\_warmup\_15\-i | Trager | 1/\(\(x\-c\)\*Sqrt\[a\*x\+b\]\) | Problem 1.15 \(i\) | \{"failed", "radicand is not a rational function over Q"\} |
| warmup\-prob\_warmup\_15\-ii | Trager | Sqrt\[a\*x\+b\]/\(x\-c\) | Problem 1.15 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_1 | Trager | \(a\*x\^2\-b\)/\(\(a\*x\^2\+b\)\*Sqrt\[a\^2\*x\^4\+b\^2\]\) | Problem 1 | \{"failed", "radicand is not a rational function over Q"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_5\-i | Trager | \(\(x\^5\+1\)\*\(4\*x\^5\-1\)\)/\(x\*\(x\^5\-a\*x\+1\)\*Sqrt\[x\^6\+x\]\) | Problem 5 \(i\) | \{"failed", "integrand is not a rational function of x and the radical over Q"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_5\-ii | Trager | \(\(x\^5\-1\)\*\(4\*x\^5\+1\)\)/\(x\*\(x\^5\-a\*x\-1\)\*Sqrt\[x\^6\-x\]\) | Problem 5 \(ii\) | \{"failed", "integrand is not a rational function of x and the radical over Q"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_7\-i | Trager | \(c\*x\^3\-d\)/\(a\*x\^3\-b\)\^\(1/3\) | Problem 7 \(i\) | \{"failed", "radicand is not a rational function over Q"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_7\-ii | Trager | \(c\*x\^3\+d\)/\(a\*x\^3\-b\)\^\(1/3\) | Problem 7 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_8\-i | Trager | \(a\*x\^4\+b\*x\^3\)\^\(1/4\) | Problem 8 \(i\) | \{"failed", "radicand is not a rational function over Q"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_8\-i | PMT | \(a\*x\^4\+b\*x\^3\)\^\(1/4\) | Problem 8 \(i\) | \{"failed", "verification failed"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_8\-ii | Trager | \(a\*x\^8\+b\*x\^7\)\^\(1/4\) | Problem 8 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_8\-ii | PMT | \(a\*x\^8\+b\*x\^7\)\^\(1/4\) | Problem 8 \(ii\) | \{"failed", "verification failed"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_8\-ii… | Trager | \(a\*x\^12\+b\*x\^11\)\^\(1/4\) | Problem 8 \(iii\) | \{"failed", "radicand is not a rational function over Q"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_8\-ii… | PMT | \(a\*x\^12\+b\*x\^11\)\^\(1/4\) | Problem 8 \(iii\) | \{"failed", "verification failed"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_9 | Trager | \(a\*x\^4\-b\)/\(x\^4\*\(2\*a\*x\^4\-b\)\^\(1/4\)\) | Problem 9 \(same as Problem 14b \(i\)\) | \{"failed", "radicand is not a rational function over Q"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_11 | Trager | \(b\*x\^4\+c\)\^4/\(x\^4\+a\*x\^2\)\^\(1/4\) | Problem 11 | \{"failed", "radicand is not a rational function over Q"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_11 | PMT | \(b\*x\^4\+c\)\^4/\(x\^4\+a\*x\^2\)\^\(1/4\) | Problem 11 | \{"failed", "partially ramified branch prime", x\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_11\-i… | PMT | \(3\*x\^4\+5\)\^4/\(x\^4\+2\*x\^2\)\^\(1/4\) | Problem 11 instantiated a=2,b=3,c=5 | \{"failed", "partially ramified branch prime", x\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_13\-i | Trager | \(a\*x\^2\+b\)/\(\(a\*x\^2\-b\)\*Sqrt\[a\*x\^3\+b\*x\^2\-b\*x\]\) | Problem 13 \(i\) | \{"failed", "radicand is not a rational function over Q"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_13\-i… | Trager | \(\(a\*x\^3\-2\*c\)\*Sqrt\[a\*x\^3\+b\*x\^2\+c\]\)/\(a\*x\^3\+c\)\^2 | Problem 13 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_13\-i… | PMT | \(\(a\*x\^3\-2\*c\)\*Sqrt\[a\*x\^3\+b\*x\^2\+c\]\)/\(a\*x\^3\+c\)\^2 | Problem 13 \(ii\) | \{"failed", "deep residue sheet undecided", c \+ a x\^3\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_13\-v… | Trager | \(\(2\*c\-a\*x\^3\)\*Sqrt\[a\*x\^3\+b\*x\^2\+c\]\)/\(\(a\*x\^3\+\(b\-3\)\*x\^2\+c\)\*\(a\*x\^3\+\(b\-2\)\*x\^2\+c\)\) | Problem 13 \(vii\) | \{"failed", "radicand is not a rational function over Q"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_14b\-… | Trager | \(x\^3\+4\*b\)/\(\(x\^3\+b\)\*\(a\*x\^4\-x\^3\-b\)\^\(1/4\)\) | Problem 14b \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_16\-i | PMT | \(\(x\^2\+1\)/\(x\*\(x\^2\-1\)\)\)\*Sqrt\[\(x\^2\-x\-1\)/\(x\^2\+x\-1\)\] | Problem 16 \(i\) | \{"needs torsion realisation \(milestone iii\)", \{\{\-1 \+ x, \{\-I, I\}\}, \{1 \+ x, \{I, \-I\}\}\}\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_16\-i… | Trager | \(\(a\*x\^2\-b\)/\(x\*\(a\*x\^2\+b\)\)\)\*Sqrt\[\(a\*x\^2\-c\*x\+b\)/\(a\*x\^2\+c\*x\+b\)\] | Problem 16 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_17\-i | PMT | \(\(2\*x\^3\-1\)/\(x\*\(x\^3\+1\)\)\)\*Sqrt\[\(x\^6\+2\*x\^3\-x\^2\+1\)/\(x\^6\+2\*x\^3\+x\^2\+1\)\] | Problem 17 \(i\) | \{"needs torsion realisation \(milestone iii\)", \{\{x, \{\-1, 1\}\}, \{1 \+ x, \{I, \-I\}\}, \{1 \- x \+ x\^2, \{I, \-I, \-I, I\}\}\}\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_17\-i… | Trager | \(\(2\*a\*x\^3\+b\)\*\(a\^2\*x\^6\-2\*a\*b\*x\^3\-a\*b\*x\^2\+b\^2\)\)/\(x\*\(a\*x\^3\-b\)\*\(a\^2\*x\^6\-2\*a\*b\*x\^3\+x\^2\+b\^2\)\)\*Sq… | Problem 17 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_17\-i… | PMT | \(\(2\*2\*x\^3\+3\)\*\(2\^2\*x\^6\-2\*2\*3\*x\^3\-2\*3\*x\^2\+3\^2\)\)/\(x\*\(2\*x\^3\-3\)\*\(2\^2\*x\^6\-2\*2\*3\*x\^3\+x\^2\+3\^2\)\)\*Sq… | Problem 17 \(ii\) instantiated a=2,b=3,c=5 | \{"needs torsion realisation \(milestone iii\)", \{\{x, \{\-1, 1\}\}, \{\-3 \+ 2 x\^3, \{\-6\*I, 6\*I, 6\*I, \-6\*I, 6\*I, \-6\*I\}\}\}\} |
| laurent\_poly\_subst\-prob\_laurent\_sub\_18\-i… | Trager | \(\(2\*a\*x\^3\+b\)/\(x\*\(a\*x\^3\+x\-b\)\)\)\*\(\(a\*x\^3\-a\*x\-b\)/\(a\*x\^3\+b\*x\-b\)\)\^\(1/3\) | Problem 18 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| rational\_subs\-prob\_rational\_sub\_1\-v | Trager | 1/\(\(x\+1\)\*\(\(a\+1\)\*x\^3\+3\*\(1\-a\)\*x\^2\+3\*\(a\+1\)\*x\-a\+1\)\^\(1/3\)\) | Problem 1 \(v\) | \{"failed", "radicand is not a rational function over Q"\} |
| rational\_subs\-prob\_rational\_sub\_1\-vi | Trager | 1/\(\(x\-1\)\*\(\(a\-1\)\*x\^3\+3\*\(a\+1\)\*x\^2\+3\*\(a\-1\)\*x\+a\+1\)\^\(1/3\)\) | Problem 1 \(vi\) | \{"failed", "radicand is not a rational function over Q"\} |
| rational\_subs\-prob\_rational\_sub\_2\-ii | PMT | \(x\^2\+2\*x\-1\)/\(x\*\(x\-1\)\*Sqrt\[2\*x\^4\-2\*x\^3\+3\*x\^2\+1\]\) | Problem 2 \(ii\) | \{"needs torsion realisation \(milestone iii\)", \{\{\-1 \+ x, \{1, \-1\}\}, \{x, \{1, \-1\}\}\}\} |
| rational\_subs\-prob\_rational\_sub\_2\-iii | Trager | \(b\*x\^2\-2\*a\*x\-a\*b\)/\(\(x\^2\+a\)\*Sqrt\[2\*x\^4\+2\*b\*x\^3\+\(2\*a\+b\^2\)\*x\^2\+a\^2\]\) | Problem 2 \(iii\) | \{"failed", "radicand is not a rational function over Q"\} |
| rational\_subs\-prob\_rational\_sub\_2\-iv | Trager | \(b\*x\^2\-2\*a\*x\-a\*b\)/\(x\*\(x\+b\)\*Sqrt\[2\*x\^4\+2\*b\*x\^3\+\(2\*a\+b\^2\)\*x\^2\+a\^2\]\) | Problem 2 \(iv\) | \{"failed", "radicand is not a rational function over Q"\} |
| rational\_subs\-prob\_rational\_sub\_2\-iv\-ins… | PMT | \(3\*x\^2\-2\*2\*x\-2\*3\)/\(x\*\(x\+3\)\*Sqrt\[2\*x\^4\+2\*3\*x\^3\+\(2\*2\+3\^2\)\*x\^2\+2\^2\]\) | Problem 2 \(iv\) instantiated a=2,b=3 | \{"needs torsion realisation \(milestone iii\)", \{\{x, \{\-1, 1\}\}, \{3 \+ x, \{\-1, 1\}\}\}\} |
| rational\_subs\-prob\_rational\_sub\_3\-ii | Trager | \(\(a\*x\^2\-2\*b\*x\-b\)\*Sqrt\[a\*\(x\^2\+x\)\*\(a\*x\^2\-a\*x\+2\*b\)\]\)/\(a\*\(b\-a\*x\)\*\(a\*x\-b\)\^2\) | Problem 3 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| rational\_subs\-prob\_rational\_sub\_3\-ii | PMT | \(\(a\*x\^2\-2\*b\*x\-b\)\*Sqrt\[a\*\(x\^2\+x\)\*\(a\*x\^2\-a\*x\+2\*b\)\]\)/\(a\*\(b\-a\*x\)\*\(a\*x\-b\)\^2\) | Problem 3 \(ii\) | \{"failed", "deep residue sheet undecided", \-b \+ a x\} |
| gunther\-prob\_gunther\_2\-ii | Trager | \(a\*x\-b\)/\(\(a\*x\+b\)\*Sqrt\[a\^2\*x\^3\+b\^2\*x\]\) | Problem 2 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| gunther\-prob\_gunther\_3\-i | Trager | x/Sqrt\[a\*x\^4\+b\*x\] | Problem 3 \(i\) | \{"failed", "radicand is not a rational function over Q"\} |
| gunther\-prob\_gunther\_3\-i | PMT | x/Sqrt\[a\*x\^4\+b\*x\] | Problem 3 \(i\) | \{"failed", "no solution within bounds", \{0\}\} |
| gunther\-prob\_gunther\_3\-i\-inst | PMT | x/Sqrt\[2\*x\^4\+3\*x\] | Problem 3 \(i\) instantiated a=2,b=3 | \{"failed", "no solution within bounds", \{0\}\} |
| gunther\-prob\_gunther\_3\-ii | Trager | \(2\*a\*x\^6\-b\)/\(\(a\*x\^6\+b\)\*Sqrt\[a\*x\^6\-c\*x\^2\+b\]\) | Problem 3 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| gunther\-prob\_gunther\_3\-iii | Trager | \(3\*a\*x\^5\-2\*b\)/\(\(a\*x\^5\-c\*x\^2\+b\)\*Sqrt\[a\*x\^5\+b\]\) | Problem 3 \(iii\) | \{"failed", "radicand is not a rational function over Q"\} |
| gunther\-prob\_gunther\_3\-iii\-inst | PMT | \(3\*2\*x\^5\-2\*3\)/\(\(2\*x\^5\-5\*x\^2\+3\)\*Sqrt\[2\*x\^5\+3\]\) | Problem 3 \(iii\) instantiated a=2,b=3,c=5 | \{"needs torsion realisation \(milestone iii\)", \{\{\-1 \+ x, \{2/Sqrt\[5\], \-2/Sqrt\[5\]\}\}, \{3 \+ 6 x \+ 4 x\^2 \+ 2 x\^3, \{\-… |
| gunther\-prob\_gunther\_3\-iv | Trager | \(a\*x\^2\+b\)/\(\(a\*x\^2\-b\)\*Sqrt\[a\*x\^3\+c\*x\^2\-b\*x\]\) | Problem 3 \(iv\) | \{"failed", "radicand is not a rational function over Q"\} |
| gunther\-prob\_gunther\_3\-v | Trager | \(4\*d\*x\^5\+b\)/\(\(a\*d\*x\^5\-c\*x\-a\*b\)\*Sqrt\[d\*x\^6\-b\*x\]\) | Problem 3 \(v\) | \{"failed", "radicand is not a rational function over Q"\} |
| gunther\-prob\_gunther\_3\-vi | Trager | \(x\^2\-1\)/\(\(x\^2\+1\)\*Sqrt\[a\*x\^4\+b\*x\^3\+c\*x\^2\+b\*x\+a\]\) | Problem 3 \(vi\) | \{"failed", "radicand is not a rational function over Q"\} |
| gunther\-prob\_gunther\_4 | Trager | \(2\*a\*x\^4\+4\*a\*x\^3\-b\*x\+b\-2\*c\)/\(\(a\*x\^4\+x\^2\+\(b\+2\)\*x\+c\+1\)\*Sqrt\[a\*x\^4\+b\*x\+c\]\) | Problem 4 | \{"failed", "radicand is not a rational function over Q"\} |
| gunther\-prob\_gunther\_5\-iii\-pp | PMT | x\^3/\(\(x\^6\+4\)\*Sqrt\[x\^6\+1\]\) | Problem 5 \(iii\) upper signs, from \(ii\) b… | \{"failed", "residue class lost: no places over the prime for a root of the residue polynomial", 4 \+ x\^6, \(\-1/1… |
| gunther\-prob\_gunther\_5\-iii\-mm | PMT | x\^3/\(\(x\^6\-4\)\*Sqrt\[x\^6\-1\]\) | Problem 5 \(iii\) lower signs, from \(ii\) b… | \{"needs torsion realisation \(milestone iii\)", \{\{\-2 \+ x\^3, \{Root\[\-1 \+ 20155392 \#1\^6 \&, 2\], Root\[\-1 \+ 20155392 \#… |
| gunther\-prob\_gunther\_5\-iv\-pm | PMT | x\^3/\(\(x\^6\+8\)\*Sqrt\[x\^6\-1\]\) | Problem 5 \(iv\) x\^6\+8 with x\^6\-1, from \(i… | \{"needs torsion realisation \(milestone iii\)", \{\{2 \+ x\^2, \{1/36\*I, \-1/36\*I, 1/36\*I, \-1/36\*I\}\}, \{4 \- 2 x\^2 \+ x\^4… |
| gunther\-prob\_gunther\_5\-iv\-mp | PMT | x\^3/\(\(x\^6\-8\)\*Sqrt\[x\^6\+1\]\) | Problem 5 \(iv\) x\^6\-8 with x\^6\+1, from \(i… | \{"needs torsion realisation \(milestone iii\)", \{\{\-2 \+ x\^2, \{1/36, \-1/36, 1/36, \-1/36\}\}, \{4 \+ 2 x\^2 \+ x\^4, \{1/64… |
| power\_expand\-power\_expand\_3 | Trager | Sqrt\[x\]/Sqrt\[x\^3\+x\^2\-x\-1\] | Problem 3 \(iii\) AXIOM example | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| generalised\_gunther\-prob\_gg\_2\-ii | Trager | \(x\-2\*a\+b\)/\(\(x\^2\+\(a\*b\^2\-2\*b\)\*x\+b\^2\-a\^2\*b\^2\)\*\(\(x\-a\)\*\(x\-b\)\)\^\(1/3\)\) | Problem gg2 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_3\-main | Trager | \(a\*x\^2\+2\*b\*x\+3\*c\)/\(\(x\^3\+a\*x\^2\+b\*x\+c\)\*\(a\*x\^2\+b\*x\+c\)\^\(1/3\)\) | Problem gg3 \(difficulty ?; elem by u=x/\(… | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_6\-i | Trager | \(x\+a\)\^2/\(\(x\+b\)\*\(2\*x\^4\+4\*\(a\+b\)\*x\^3\+6\*\(a\^2\+b\^2\)\*x\^2\+4\*\(a\^3\+b\^3\)\*x\+a\^4\+b\^4\)\^\(3/4\)\) | Problem gg6 \(i\) | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_6\-ii | Trager | \(x\-a\)\^2/\(\(x\-b\)\*\(2\*x\^4\-4\*\(a\+b\)\*x\^3\+6\*\(a\^2\+b\^2\)\*x\^2\-4\*\(a\^3\+b\^3\)\*x\+a\^4\+b\^4\)\^\(3/4\)\) | Problem gg6 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_9\-i | Trager | x/\(\(2\*a\*x\^3\+b\)\*\(a\*x\^3\+b\)\^\(1/3\)\) | Problem gg9 \(i\) | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_9\-i | PMT | x/\(\(2\*a\*x\^3\+b\)\*\(a\*x\^3\+b\)\^\(1/3\)\) | Problem gg9 \(i\) | \{"failed", "no solution within bounds", \{0\}\} |
| generalised\_gunther\-prob\_gg\_9\-ii | Trager | x\^3/\(\(2\*a\*x\^3\+b\)\*\(a\*x\^3\+b\)\^\(1/3\)\) | Problem gg9 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_9\-ii | PMT | x\^3/\(\(2\*a\*x\^3\+b\)\*\(a\*x\^3\+b\)\^\(1/3\)\) | Problem gg9 \(ii\) | \{"failed", "no solution within bounds", \{0\}\} |
| generalised\_gunther\-prob\_gg\_9\-iii\_n1 | Trager | \(2\*a\*x\^3\+b\)\*x/\(\(a\^2\*m\*x\^6\+a\*b\*m\*x\^3\+b\^2\)\*\(a\*x\^3\+b\)\^\(1/3\)\) | Problem gg9 \(iii\) general form with n=1 | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_9\-iii\_n2 | Trager | \(3\*a\*x\^3\+2\*b\)\*x\^3/\(\(a\^2\*m\*x\^9\+a\*b\*m\*x\^6\+b\^2\)\*\(a\*x\^3\+b\)\^\(1/3\)\) | Problem gg9 \(iii\) general form with n=2 | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_9\-iii\_n3 | Trager | \(4\*a\*x\^3\+3\*b\)\*x\^5/\(\(a\^2\*m\*x\^12\+a\*b\*m\*x\^9\+b\^2\)\*\(a\*x\^3\+b\)\^\(1/3\)\) | Problem gg9 \(iii\) general form with n=3 | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_10\-n2 | Trager | 1/\(\(2\*b\*c/a\*x\^2\+c\)\*Sqrt\[a\+b\*x\^2\]\) | Problem gg10 Albert Rich form x\^\(n\-2\)/\(\(… | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_10\-n3 | Trager | x/\(\(2\*b\*c/a\*x\^3\+c\)\*\(a\+b\*x\^3\)\^\(1/3\)\) | Problem gg10 Albert Rich form x\^\(n\-2\)/\(\(… | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_10\-n3 | PMT | x/\(\(2\*b\*c/a\*x\^3\+c\)\*\(a\+b\*x\^3\)\^\(1/3\)\) | Problem gg10 Albert Rich form x\^\(n\-2\)/\(\(… | \{"failed", "no solution within bounds", \{0\}\} |
| generalised\_gunther\-prob\_gg\_10\-n4 | Trager | x\^2/\(\(2\*b\*c/a\*x\^4\+c\)\*\(a\+b\*x\^4\)\^\(1/4\)\) | Problem gg10 Albert Rich form x\^\(n\-2\)/\(\(… | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_11\-ii\_plus | PMT | 1/\(\(x\+1\)\*\(x\^2\+2\*x\^3\)\^\(1/4\)\) | Problem gg11 \(ii\) signs \+ \+ | \{"failed", "partially ramified branch prime", x\} |
| generalised\_gunther\-prob\_gg\_11\-iii\_plus | Trager | 1/\(\(a\*x\+b\)\*\(a\^2\*x\^3\+b\^2\*x\)\^\(1/4\)\) | Problem gg11 \(iii\) sign \+ | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_11\-iii\_minus | Trager | 1/\(\(a\*x\-b\)\*\(a\^2\*x\^3\+b\^2\*x\)\^\(1/4\)\) | Problem gg11 \(iii\) sign \- | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_11\-iv\_plus | Trager | 1/\(\(a\*x\+b\)\*\(b\*x\^2\+2\*a\*x\^3\)\^\(1/4\)\) | Problem gg11 \(iv\) signs \+ \+ | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_11\-iv\_plus | PMT | 1/\(\(a\*x\+b\)\*\(b\*x\^2\+2\*a\*x\^3\)\^\(1/4\)\) | Problem gg11 \(iv\) signs \+ \+ | \{"failed", "partially ramified branch prime", x\} |
| generalised\_gunther\-prob\_gg\_11\-iv\_plus\-i… | PMT | 1/\(\(2\*x\+3\)\*\(3\*x\^2\+2\*2\*x\^3\)\^\(1/4\)\) | Problem gg11 \(iv\) signs \+ \+ instantiated… | \{"failed", "partially ramified branch prime", x\} |
| generalised\_gunther\-prob\_gg\_11\-iv\_minus | Trager | 1/\(\(a\*x\-b\)\*\(b\*x\^2\-2\*a\*x\^3\)\^\(1/4\)\) | Problem gg11 \(iv\) signs \- \- | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_11\-iv\_minus | PMT | 1/\(\(a\*x\-b\)\*\(b\*x\^2\-2\*a\*x\^3\)\^\(1/4\)\) | Problem gg11 \(iv\) signs \- \- | \{"failed", "partially ramified branch prime", x\} |
| generalised\_gunther\-prob\_gg\_11\-iv\_minus\-… | PMT | 1/\(\(2\*x\-3\)\*\(3\*x\^2\-2\*2\*x\^3\)\^\(1/4\)\) | Problem gg11 \(iv\) signs \- \- instantiated… | \{"failed", "partially ramified branch prime", x\} |
| generalised\_gunther\-prob\_gg\_13\-i | Trager | \(a\*c\*x\^4\+2\*b\*c\*x\^2\+b\*d\)\*Sqrt\[\(a\*x\^2\+b\)/\(c\*x\^2\+d\)\]/\(c\^2\*x\^8\+2\*c\*d\*x\^6\+\(d\^2\-a\^2\)\*x\^4\-2\*a\*b\*x… | Problem gg13 \(i\) rational radicand | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_13\-i | PMT | \(a\*c\*x\^4\+2\*b\*c\*x\^2\+b\*d\)\*Sqrt\[\(a\*x\^2\+b\)/\(c\*x\^2\+d\)\]/\(c\^2\*x\^8\+2\*c\*d\*x\^6\+\(d\^2\-a\^2\)\*x\^4\-2\*a\*b\*x… | Problem gg13 \(i\) rational radicand | \{"needs torsion realisation \(milestone iii\)", \{\{b \+ a x\^2 \- d x\^2 \- c x\^4, \{\-\(2 a\^3 \+ 4 b c Sqrt\[4 b c \+ \(a \- … |
| generalised\_gunther\-prob\_gg\_13\-ii | Trager | \(2\*a\*c\*x\^4\+\(a\*d\+3\*b\*c\)\*x\^2\+2\*b\*d\)\*\(\(a\*x\^2\+b\)/\(c\*x\^2\+d\)\)\^\(3/4\)/\(\(a\*x\^2\+b\)\*\(c\*x\^6\+d\*x\^4\-a\*x\^… | Problem gg13 \(ii\) rational radicand | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_13\-iii | Trager | \(2\*a\*c\*x\^3\+a\*d\*x\^2\+4\*b\*c\*x\+3\*b\*d\)/\(\(c\*x\+d\)\*\(c\*x\^4\+d\*x\^3\-a\*x\^2\-b\)\*\(\(a\*x\^2\+b\)/\(c\*x\+d\)\)\^\(1/3\)… | Problem gg13 \(iii\) rational radicand | \{"failed", "radicand is not a rational function over Q"\} |
| generalised\_gunther\-prob\_gg\_13\-iv | Trager | \(a\*c\*x\^6\+2\*b\*c\*x\^3\+b\*d\)\*\(\(a\*x\^3\+b\)/\(c\*x\^3\+d\)\)\^\(2/3\)/\(\(a\*x\^3\+b\)\*\(c\*x\^6\+\(d\-a\)\*x\^3\-b\)\) | Problem gg13 \(iv\) rational radicand | \{"failed", "radicand is not a rational function over Q"\} |
| chebychev\-intro\-1 | Trager | x\^\(\-5/6\)\*Sqrt\[x\^\(1/3\)\-1\] | Chapter intro worked example \(type III\) | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_1\-i | Trager | x\^\(1/3\)\*\(a\*Sqrt\[x\]\+b\)\^3 | Problem chebychev\_1 \(i\) | \{"failed", "integrand is not a rational function of x and the radical over Q"\} |
| chebychev\-chebychev\_1\-ii | Trager | x\^\(3/2\)\*\(a\*x\^\(5/6\)\+b\)\^\(1/3\) | Problem chebychev\_1 \(ii\) | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_1\-ii | PMT | x\^\(3/2\)\*\(a\*x\^\(5/6\)\+b\)\^\(1/3\) | Problem chebychev\_1 \(ii\) | \{"failed", "verification failed"\} |
| chebychev\-chebychev\_1\-ii\-inst | Trager | x\^\(3/2\)\*\(2\*x\^\(5/6\)\+3\)\^\(1/3\) | Problem chebychev\_1 \(ii\) instantiated a=… | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_1\-iii | Trager | x\^\(4/3\)/\(a\*Sqrt\[x\]\+b\)\^\(2/3\) | Problem chebychev\_1 \(iii\) | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_1\-iii\-inst | Trager | x\^\(4/3\)/\(2\*Sqrt\[x\]\+3\)\^\(2/3\) | Problem chebychev\_1 \(iii\) instantiated a… | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_1\-iv | Trager | \(a\*x\^\(3/4\)\+b\)\^\(2/3\)/x\^3 | Problem chebychev\_1 \(iv\) | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_1\-iv\-inst | Trager | \(2\*x\^\(3/4\)\+3\)\^\(2/3\)/x\^3 | Problem chebychev\_1 \(iv\) instantiated a=… | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_3\-i | Trager | \(a\*x\^\(1/5\)\+b\)\^\(3/4\) | Problem chebychev\_3 \(i\) | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_3\-i\-inst | Trager | \(2\*x\^\(1/5\)\+3\)\^\(3/4\) | Problem chebychev\_3 \(i\) instantiated a=2… | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_3\-iii | Trager | x\^\(1/3\)\*\(a\*x\^2\+b\)\^\(1/3\) | Problem chebychev\_3 \(iii\) | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| chebychev\-chebychev\_3\-iii | PMT | x\^\(1/3\)\*\(a\*x\^2\+b\)\^\(1/3\) | Problem chebychev\_3 \(iii\) | \{"failed", "internal: solution does not verify"\} |
| chebychev\-chebychev\_3\-iii\-inst | Trager | x\^\(1/3\)\*\(2\*x\^2\+3\)\^\(1/3\) | Problem chebychev\_3 \(iii\) instantiated a… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| chebychev\-chebychev\_3\-v | Trager | 1/\(a\*x\^3\+b\)\^\(1/3\) | Problem chebychev\_3 \(v\) | \{"failed", "radicand is not a rational function over Q"\} |
| chebychev\-chebychev\_3\-v | PMT | 1/\(a\*x\^3\+b\)\^\(1/3\) | Problem chebychev\_3 \(v\) | \{"failed", "internal: solution does not verify"\} |
| chebychev\-chebychev\_3\-vii | Trager | x\^4/\(a\*x\^4\+b\)\^\(5/4\) | Problem chebychev\_3 \(vii\) \(hint: q\+\(p\+1\)… | \{"failed", "radicand is not a rational function over Q"\} |
| chebychev\-chebychev\_3\-vii | PMT | x\^4/\(a\*x\^4\+b\)\^\(5/4\) | Problem chebychev\_3 \(vii\) \(hint: q\+\(p\+1\)… | \{"failed", "no solution within bounds", \{4\}\} |
| chebychev\-chebychev\_3\-viii | Trager | x\^\(4/5\)\*\(a/x\^\(4/5\)\+b\)\^\(5/4\) | Problem chebychev\_3 \(viii\) \(expect by bo… | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_3\-viii\-inst | Trager | x\^\(4/5\)\*\(2/x\^\(4/5\)\+3\)\^\(5/4\) | Problem chebychev\_3 \(viii\) \(expect by bo… | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_3\-x | Trager | x\^\(3/2\)/Sqrt\[a/Sqrt\[x\]\+b\] | Problem chebychev\_3 \(x\) \(expect by book … | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_3\-x\-inst | Trager | x\^\(3/2\)/Sqrt\[2/Sqrt\[x\]\+3\] | Problem chebychev\_3 \(x\) \(expect by book … | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_3b\-v | Trager | Sqrt\[Tan\[x\]\] | Problem chebychev\_3b \(v\) | \{"failed", "transcendental function in the integrand"\} |
| chebychev\-chebychev\_3b\-vi | Trager | Sqrt\[\(\(Exp\[x\]\-Exp\[\-x\]\)/\(Exp\[x\]\+Exp\[\-x\]\)\)\] | Problem chebychev\_3b \(vi\) tanh via Exp | \{"failed", "transcendental function in the integrand"\} |
| chebychev\-chebychev\_3b\-vii | Trager | Sqrt\[Cos\[x\]/Sin\[x\]\] | Problem chebychev\_3b \(vii\) cot as Cos/Si… | \{"failed", "transcendental function in the integrand"\} |
| chebychev\-chebychev\_3b\-viii | Trager | Sqrt\[\(\(Exp\[x\]\+Exp\[\-x\]\)/\(Exp\[x\]\-Exp\[\-x\]\)\)\] | Problem chebychev\_3b \(viii\) coth via Exp | \{"failed", "transcendental function in the integrand"\} |
| chebychev\-chebychev\_3c\-iii\_n2 | Trager | 1/\(Tan\[x\]\)\^\(1/2\) | Problem chebychev\_3c \(iii\) 1/tan\^\(1/n\) n… | \{"failed", "transcendental function in the integrand"\} |
| chebychev\-chebychev\_3c\-vi\_n2 | Trager | \(Tan\[x\]\)\^\(1/2\) | Problem chebychev\_3c \(vi\) 1/cot\^\(1/n\)=ta… | \{"failed", "transcendental function in the integrand"\} |
| chebychev\-chebychev\_3c\-iii\_n3 | Trager | 1/\(Tan\[x\]\)\^\(1/3\) | Problem chebychev\_3c \(iii\) 1/tan\^\(1/n\) n… | \{"failed", "transcendental function in the integrand"\} |
| chebychev\-chebychev\_3c\-vi\_n3 | Trager | \(Tan\[x\]\)\^\(1/3\) | Problem chebychev\_3c \(vi\) 1/cot\^\(1/n\)=ta… | \{"failed", "transcendental function in the integrand"\} |
| chebychev\-chebychev\_2\-i | Trager | \(x\^2\-2\*x\-1\)\*\(2\*Sqrt\[x\-1\]/Sqrt\[x\^2\+1\]\-1\)\^\(2/3\)/\(\(x\^2\+1\)\^2\*\(Sqrt\[x\-1\]/Sqrt\[x\^2\+1\]\)\^\(2/3\)\) | Problem chebychev\_2 \(i\) \(solved\) | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_2\-i | PMT | \(x\^2\-2\*x\-1\)\*\(2\*Sqrt\[x\-1\]/Sqrt\[x\^2\+1\]\-1\)\^\(2/3\)/\(\(x\^2\+1\)\^2\*\(Sqrt\[x\-1\]/Sqrt\[x\^2\+1\]\)\^\(2/3\)\) | Problem chebychev\_2 \(i\) \(solved\) | \{"failed", "tower construction failed", \(\(\-1 \- 2 x \+ x\^2\) \(\-1 \+ \(2 Sqrt\[\-1 \+ x\]\)/Sqrt\[1 \+ x\^2\]\)\^\(2/3\)\)/\(\(Sqrt\[… |
| chebychev\-chebychev\_2\-ii | Trager | \(x\^3\-x\^2\-x\-1\)\*\(\(Sqrt\[x\-1\]/\(Sqrt\[x\+1\]\*Sqrt\[x\^2\+1\]\)\)\-1\)\^\(2/3\)/\(Sqrt\[x\-1\]\*\(x\+1\)\^\(3/2\)\*\(x\^2\+1\)… | Problem chebychev\_2 \(ii\) \(hint substitut… | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_2\-ii | PMT | \(x\^3\-x\^2\-x\-1\)\*\(\(Sqrt\[x\-1\]/\(Sqrt\[x\+1\]\*Sqrt\[x\^2\+1\]\)\)\-1\)\^\(2/3\)/\(Sqrt\[x\-1\]\*\(x\+1\)\^\(3/2\)\*\(x\^2\+1\)… | Problem chebychev\_2 \(ii\) \(hint substitut… | \{"failed", "tower construction failed", \(\(\-1 \- x \- x\^2 \+ x\^3\) \(\-1 \+ Sqrt\[\-1 \+ x\]/\(Sqrt\[1 \+ x\] Sqrt\[1 \+ x\^2\]\)\)\^… |
| chebychev\-chebychev\_5\-i | Trager | Sqrt\[a\*x\^2\+b\]/Sqrt\[c/x\^2\+d\] | Problem chebychev\_5 \(i\) type III conject… | \{"failed", "radicand is not a rational function over Q"\} |
| chebychev\-chebychev\_5\-i | PMT | Sqrt\[a\*x\^2\+b\]/Sqrt\[c/x\^2\+d\] | Problem chebychev\_5 \(i\) type III conject… | \{"failed", "tower construction failed", Sqrt\[b \+ a x\^2\]/Sqrt\[d \+ c/x\^2\]\} |
| chebychev\-chebychev\_5\-i\-inst | Trager | Sqrt\[2\*x\^2\+3\]/Sqrt\[5/x\^2\+7\] | Problem chebychev\_5 \(i\) type III conject… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| chebychev\-chebychev\_5\-i\-inst | PMT | Sqrt\[2\*x\^2\+3\]/Sqrt\[5/x\^2\+7\] | Problem chebychev\_5 \(i\) type III conject… | \{"failed", "tower construction failed", Sqrt\[3 \+ 2 x\^2\]/Sqrt\[7 \+ 5/x\^2\]\} |
| chebychev\-chebychev\_5\-ii | Trager | Sqrt\[a/x\^2\+b\]/Sqrt\[c\*x\^2\+d\] | Problem chebychev\_5 \(ii\) type III conjec… | \{"failed", "radicand is not a rational function over Q"\} |
| chebychev\-chebychev\_5\-ii | PMT | Sqrt\[a/x\^2\+b\]/Sqrt\[c\*x\^2\+d\] | Problem chebychev\_5 \(ii\) type III conjec… | \{"failed", "tower construction failed", Sqrt\[b \+ a/x\^2\]/Sqrt\[d \+ c x\^2\]\} |
| chebychev\-chebychev\_5\-ii\-inst | Trager | Sqrt\[2/x\^2\+3\]/Sqrt\[5\*x\^2\+7\] | Problem chebychev\_5 \(ii\) type III conjec… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| chebychev\-chebychev\_5\-ii\-inst | PMT | Sqrt\[2/x\^2\+3\]/Sqrt\[5\*x\^2\+7\] | Problem chebychev\_5 \(ii\) type III conjec… | \{"failed", "tower construction failed", Sqrt\[3 \+ 2/x\^2\]/Sqrt\[7 \+ 5 x\^2\]\} |
| chebychev\-chebychev\_5\-iii | Trager | Sqrt\[a\*x\^2\+b\]\*Sqrt\[c/x\^2\+d\] | Problem chebychev\_5 \(iii\) type III conje… | \{"failed", "radicand is not a rational function over Q"\} |
| chebychev\-chebychev\_5\-iii | PMT | Sqrt\[a\*x\^2\+b\]\*Sqrt\[c/x\^2\+d\] | Problem chebychev\_5 \(iii\) type III conje… | \{"failed", "tower construction failed", Sqrt\[d \+ c/x\^2\] Sqrt\[b \+ a x\^2\]\} |
| chebychev\-chebychev\_5\-iii\-inst | Trager | Sqrt\[2\*x\^2\+3\]\*Sqrt\[5/x\^2\+7\] | Problem chebychev\_5 \(iii\) type III conje… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| chebychev\-chebychev\_5\-iii\-inst | PMT | Sqrt\[2\*x\^2\+3\]\*Sqrt\[5/x\^2\+7\] | Problem chebychev\_5 \(iii\) type III conje… | \{"failed", "tower construction failed", Sqrt\[7 \+ 5/x\^2\] Sqrt\[3 \+ 2 x\^2\]\} |
| chebychev\-chebychev\_5\-iv | Trager | 1/\(Sqrt\[a\*x\^2\+b\]\*Sqrt\[c/x\^2\+d\]\) | Problem chebychev\_5 \(iv\) type III conjec… | \{"failed", "radicand is not a rational function over Q"\} |
| chebychev\-chebychev\_5\-iv | PMT | 1/\(Sqrt\[a\*x\^2\+b\]\*Sqrt\[c/x\^2\+d\]\) | Problem chebychev\_5 \(iv\) type III conjec… | \{"failed", "tower construction failed", 1/\(Sqrt\[d \+ c/x\^2\] Sqrt\[b \+ a x\^2\]\)\} |
| chebychev\-chebychev\_5\-iv\-inst | Trager | 1/\(Sqrt\[2\*x\^2\+3\]\*Sqrt\[5/x\^2\+7\]\) | Problem chebychev\_5 \(iv\) type III conjec… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| chebychev\-chebychev\_5\-iv\-inst | PMT | 1/\(Sqrt\[2\*x\^2\+3\]\*Sqrt\[5/x\^2\+7\]\) | Problem chebychev\_5 \(iv\) type III conjec… | \{"failed", "tower construction failed", 1/\(Sqrt\[7 \+ 5/x\^2\] Sqrt\[3 \+ 2 x\^2\]\)\} |
| chebychev\-chebychev\_5\-v | Trager | \(a\*Sqrt\[x\]\+b\)\^\(2/3\)/\(c\*Sqrt\[x\]\+d\)\^\(2/3\) | Problem chebychev\_5 \(v\) type IV conjectu… | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_5\-v\-inst | Trager | \(2\*Sqrt\[x\]\+3\)\^\(2/3\)/\(5\*Sqrt\[x\]\+7\)\^\(2/3\) | Problem chebychev\_5 \(v\) type IV conjectu… | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_5\-vi | Trager | \(a/x\^\(1/3\)\+b\)\^\(1/4\)/\(c/x\^\(1/3\)\+d\)\^\(1/4\) | Problem chebychev\_5 \(vi\) type IV conject… | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_5\-vi | PMT | \(a/x\^\(1/3\)\+b\)\^\(1/4\)/\(c/x\^\(1/3\)\+d\)\^\(1/4\) | Problem chebychev\_5 \(vi\) type IV conject… | \{"failed", "tower construction failed", \(b \+ a/x\^\(1/3\)\)\^\(1/4\)/\(d \+ c/x\^\(1/3\)\)\^\(1/4\)\} |
| chebychev\-chebychev\_5\-vi\-inst | Trager | \(2/x\^\(1/3\)\+3\)\^\(1/4\)/\(5/x\^\(1/3\)\+7\)\^\(1/4\) | Problem chebychev\_5 \(vi\) type IV conject… | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_5\-vi\-inst | PMT | \(2/x\^\(1/3\)\+3\)\^\(1/4\)/\(5/x\^\(1/3\)\+7\)\^\(1/4\) | Problem chebychev\_5 \(vi\) type IV conject… | \{"failed", "tower construction failed", \(3 \+ 2/x\^\(1/3\)\)\^\(1/4\)/\(7 \+ 5/x\^\(1/3\)\)\^\(1/4\)\} |
| chebychev\-chebychev\_5\-vii | Trager | \(a/x\+b\)\^\(1/3\)\*\(c/x\+d\)\^\(2/3\) | Problem chebychev\_5 \(vii\) type V conject… | \{"failed", "radicand is not a rational function over Q"\} |
| chebychev\-chebychev\_5\-vii | PMT | \(a/x\+b\)\^\(1/3\)\*\(c/x\+d\)\^\(2/3\) | Problem chebychev\_5 \(vii\) type V conject… | \{"failed", "tower construction failed", \(b \+ a/x\)\^\(1/3\) \(d \+ c/x\)\^\(2/3\)\} |
| chebychev\-chebychev\_5\-vii\-inst | Trager | \(2/x\+3\)\^\(1/3\)\*\(5/x\+7\)\^\(2/3\) | Problem chebychev\_5 \(vii\) type V conject… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| chebychev\-chebychev\_5\-vii\-inst | PMT | \(2/x\+3\)\^\(1/3\)\*\(5/x\+7\)\^\(2/3\) | Problem chebychev\_5 \(vii\) type V conject… | \{"failed", "tower construction failed", \(3 \+ 2/x\)\^\(1/3\) \(7 \+ 5/x\)\^\(2/3\)\} |
| chebychev\-chebychev\_5\-viii | Trager | \(a\*x\^\(1/3\)\+b\)\^\(5/3\)/\(c\*x\^\(1/3\)\+d\)\^\(2/3\) | Problem chebychev\_5 \(viii\) type V conjec… | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_5\-viii\-inst | Trager | \(2\*x\^\(1/3\)\+3\)\^\(5/3\)/\(5\*x\^\(1/3\)\+7\)\^\(2/3\) | Problem chebychev\_5 \(viii\) type V conjec… | \{"failed", "nested radicals"\} |
| chebychev\-chebychev\_5\-final\_special\_n2 | Trager | \(a\*x\^2\+b\)\^\(1/2\)/\(b\^2\-a\^2\*x\^4\)\^\(1/2\) | Problem chebychev\_5 closing remark, n=2,… | \{"failed", "radicand is not a rational function over Q"\} |
| chebychev\-chebychev\_5\-final\_special\_n2 | PMT | \(a\*x\^2\+b\)\^\(1/2\)/\(b\^2\-a\^2\*x\^4\)\^\(1/2\) | Problem chebychev\_5 closing remark, n=2,… | \{"failed", "tower construction failed", Sqrt\[b \+ a x\^2\]/Sqrt\[b\^2 \- a\^2 x\^4\]\} |
| chebychev\-chebychev\_5\-final\_special\_n2\-i… | Trager | \(2\*x\^2\+3\)\^\(1/2\)/\(3\^2\-2\^2\*x\^4\)\^\(1/2\) | Problem chebychev\_5 closing remark, n=2,… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| chebychev\-chebychev\_5\-final\_special\_n2\-i… | PMT | \(2\*x\^2\+3\)\^\(1/2\)/\(3\^2\-2\^2\*x\^4\)\^\(1/2\) | Problem chebychev\_5 closing remark, n=2,… | \{"failed", "tower construction failed", Sqrt\[3 \+ 2 x\^2\]/Sqrt\[9 \- 4 x\^4\]\} |
| chebychev\-chebychev\_5\-final\_special\_n3 | Trager | \(a\*x\^3\+b\)\^\(1/3\)/\(b\^2\-a\^2\*x\^6\)\^\(1/3\) | Problem chebychev\_5 closing remark, n=3,… | \{"failed", "radicand is not a rational function over Q"\} |
| chebychev\-chebychev\_5\-final\_special\_n3 | PMT | \(a\*x\^3\+b\)\^\(1/3\)/\(b\^2\-a\^2\*x\^6\)\^\(1/3\) | Problem chebychev\_5 closing remark, n=3,… | \{"failed", "tower construction failed", \(b \+ a x\^3\)\^\(1/3\)/\(b\^2 \- a\^2 x\^6\)\^\(1/3\)\} |
| chebychev\-chebychev\_5\-final\_special\_n3\-i… | Trager | \(2\*x\^3\+3\)\^\(1/3\)/\(3\^2\-2\^2\*x\^6\)\^\(1/3\) | Problem chebychev\_5 closing remark, n=3,… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| chebychev\-chebychev\_5\-final\_special\_n3\-i… | PMT | \(2\*x\^3\+3\)\^\(1/3\)/\(3\^2\-2\^2\*x\^6\)\^\(1/3\) | Problem chebychev\_5 closing remark, n=3,… | \{"failed", "tower construction failed", \(3 \+ 2 x\^3\)\^\(1/3\)/\(9 \- 4 x\^6\)\^\(1/3\)\} |
| inverse\_function\-intro\-1 | Trager | Sqrt\[1\-3/x\+Sqrt\[9\*x\^2\-6\*x\+9\]/x\] | Chapter intro worked example \(solved\) | \{"failed", "nested radicals"\} |
| inverse\_function\-intro\-1 | PMT | Sqrt\[1\-3/x\+Sqrt\[9\*x\^2\-6\*x\+9\]/x\] | Chapter intro worked example \(solved\) | \{"failed", "tower construction failed", Sqrt\[1 \- 3/x \+ Sqrt\[9 \- 6 x \+ 9 x\^2\]/x\]\} |
| inverse\_function\-inverse\_function\_1\-main | Trager | Sqrt\[a\-b\-a/x\-b/x\+Sqrt\[\(a\+b\)\^2\*x\^2\-2\*\(a\^2\-b\^2\)\*x\+\(a\+b\)\^2\]/x\] | Problem inverse\_function\_1 | \{"failed", "nested radicals"\} |
| inverse\_function\-inverse\_function\_1\-main | PMT | Sqrt\[a\-b\-a/x\-b/x\+Sqrt\[\(a\+b\)\^2\*x\^2\-2\*\(a\^2\-b\^2\)\*x\+\(a\+b\)\^2\]/x\] | Problem inverse\_function\_1 | \{"failed", "tower construction failed", Sqrt\[a \- b \- a/x \- b/x \+ Sqrt\[\(a \+ b\)\^2 \- 2 \(a\^2 \- b\^2\) x \+ \(a \+ b\)\^2 … |
| inverse\_function\-inverse\_function\_1\-main… | Trager | Sqrt\[2\-3\-2/x\-3/x\+Sqrt\[\(2\+3\)\^2\*x\^2\-2\*\(2\^2\-3\^2\)\*x\+\(2\+3\)\^2\]/x\] | Problem inverse\_function\_1 instantiated … | \{"failed", "nested radicals"\} |
| inverse\_function\-inverse\_function\_1\-main… | PMT | Sqrt\[2\-3\-2/x\-3/x\+Sqrt\[\(2\+3\)\^2\*x\^2\-2\*\(2\^2\-3\^2\)\*x\+\(2\+3\)\^2\]/x\] | Problem inverse\_function\_1 instantiated … | \{"failed", "tower construction failed", Sqrt\[\-1 \- 5/x \+ Sqrt\[25 \+ 10 x \+ 25 x\^2\]/x\]\} |
| inverse\_function\-inverse\_function\_2\-i | Trager | Sqrt\[\(a\-b\)\*x\+Sqrt\[a\+b\]\*Sqrt\[x\]\*Sqrt\[\(a\+b\)\*x\+4\*a\*b\]\]/Sqrt\[x\+a\+b\] | Problem inverse\_function\_2 \(i\) | \{"failed", "algebraic constant in the integrand \(coefficients must be rational\)"\} |
| inverse\_function\-inverse\_function\_2\-i | PMT | Sqrt\[\(a\-b\)\*x\+Sqrt\[a\+b\]\*Sqrt\[x\]\*Sqrt\[\(a\+b\)\*x\+4\*a\*b\]\]/Sqrt\[x\+a\+b\] | Problem inverse\_function\_2 \(i\) | \{"failed", "tower construction failed", Sqrt\[\(a \- b\) x \+ Sqrt\[a \+ b\] Sqrt\[x\] Sqrt\[4 a b \+ \(a \+ b\) x\]\]/Sqrt\[a \+… |
| inverse\_function\-inverse\_function\_2\-i\-in… | Trager | Sqrt\[\(2\-3\)\*x\+Sqrt\[2\+3\]\*Sqrt\[x\]\*Sqrt\[\(2\+3\)\*x\+4\*2\*3\]\]/Sqrt\[x\+2\+3\] | Problem inverse\_function\_2 \(i\) instantia… | \{"failed", "algebraic constant in the integrand \(coefficients must be rational\)"\} |
| inverse\_function\-inverse\_function\_2\-i\-in… | PMT | Sqrt\[\(2\-3\)\*x\+Sqrt\[2\+3\]\*Sqrt\[x\]\*Sqrt\[\(2\+3\)\*x\+4\*2\*3\]\]/Sqrt\[x\+2\+3\] | Problem inverse\_function\_2 \(i\) instantia… | \{"failed", "tower construction failed", Sqrt\[\-x \+ Sqrt\[5\] Sqrt\[x\] Sqrt\[24 \+ 5 x\]\]/Sqrt\[5 \+ x\]\} |
| inverse\_function\-inverse\_function\_2\-ii | Trager | \(\(a\-b\)\*x\-Sqrt\[x\]\*Sqrt\[\(a\+b\)\^2\*x\+4\*a\*b\]\)\^\(1/3\)/\(x\+1\)\^\(1/3\) | Problem inverse\_function\_2 \(ii\) | \{"failed", "nested radicals"\} |
| inverse\_function\-inverse\_function\_2\-ii | PMT | \(\(a\-b\)\*x\-Sqrt\[x\]\*Sqrt\[\(a\+b\)\^2\*x\+4\*a\*b\]\)\^\(1/3\)/\(x\+1\)\^\(1/3\) | Problem inverse\_function\_2 \(ii\) | \{"failed", "tower construction failed", \(\(a \- b\) x \- Sqrt\[x\] Sqrt\[4 a b \+ \(a \+ b\)\^2 x\]\)\^\(1/3\)/\(1 \+ x\)\^\(1/3\)\} |
| inverse\_function\-inverse\_function\_2\-ii\-i… | Trager | \(\(2\-3\)\*x\-Sqrt\[x\]\*Sqrt\[\(2\+3\)\^2\*x\+4\*2\*3\]\)\^\(1/3\)/\(x\+1\)\^\(1/3\) | Problem inverse\_function\_2 \(ii\) instanti… | \{"failed", "nested radicals"\} |
| inverse\_function\-inverse\_function\_2\-ii\-i… | PMT | \(\(2\-3\)\*x\-Sqrt\[x\]\*Sqrt\[\(2\+3\)\^2\*x\+4\*2\*3\]\)\^\(1/3\)/\(x\+1\)\^\(1/3\) | Problem inverse\_function\_2 \(ii\) instanti… | \{"failed", "tower construction failed", \(\-x \- Sqrt\[x\] Sqrt\[24 \+ 25 x\]\)\^\(1/3\)/\(1 \+ x\)\^\(1/3\)\} |
| inverse\_function\-inverse\_function\_2\-iii | Trager | \(\(Sqrt\[x\^3\+81\*x\^2\]\+9\*x\)\^\(2/3\)\-x\)\^\(1/3\)/\(Sqrt\[x\^3\+81\*x\^2\]\+9\*x\)\^\(1/9\) | Problem inverse\_function\_2 \(iii\) | \{"failed", "nested radicals"\} |
| inverse\_function\-inverse\_function\_2\-iii | PMT | \(\(Sqrt\[x\^3\+81\*x\^2\]\+9\*x\)\^\(2/3\)\-x\)\^\(1/3\)/\(Sqrt\[x\^3\+81\*x\^2\]\+9\*x\)\^\(1/9\) | Problem inverse\_function\_2 \(iii\) | \{"failed", "tower construction failed", \(\-x \+ \(9 x \+ Sqrt\[81 x\^2 \+ x\^3\]\)\^\(2/3\)\)\^\(1/3\)/\(9 x \+ Sqrt\[81 x\^2 \+ x\^3… |
| inverse\_function\-inverse\_function\_2\-iv | Trager | \(3\^\(1/3\)\*x\+\(Sqrt\[3\]\*Sqrt\[27\-x\^3\]\+9\)\^\(2/3\)\)/\(Sqrt\[3\]\*Sqrt\[27\-x\^3\]\+9\)\^\(1/3\) | Problem inverse\_function\_2 \(iv\) | \{"failed", "algebraic constant in the integrand \(coefficients must be rational\)"\} |
| inverse\_function\-inverse\_function\_2\-iv | PMT | \(3\^\(1/3\)\*x\+\(Sqrt\[3\]\*Sqrt\[27\-x\^3\]\+9\)\^\(2/3\)\)/\(Sqrt\[3\]\*Sqrt\[27\-x\^3\]\+9\)\^\(1/3\) | Problem inverse\_function\_2 \(iv\) | \{"failed", "tower construction failed", \(3\^\(1/3\) x \+ \(9 \+ Sqrt\[3\] Sqrt\[27 \- x\^3\]\)\^\(2/3\)\)/\(9 \+ Sqrt\[3\] Sqrt\[27 … |
| inverse\_function\-inverse\_function\_2\-v | Trager | \(2\^\(1/3\)\*\(Sqrt\[x\^2\+4\]\+x\)\^\(2/3\)\-2\)\^\(1/3\)/\(Sqrt\[x\^2\+4\]\+x\)\^\(1/9\) | Problem inverse\_function\_2 \(v\) | \{"failed", "algebraic constant in the integrand \(coefficients must be rational\)"\} |
| inverse\_function\-inverse\_function\_3\-iii | Trager | ArcTan\[x\^2\] | Problem inverse\_function\_3 \(iii\) BOOK CL… | \{"failed", "transcendental function in the integrand"\} |
| inverse\_function\-inverse\_function\_3\-vi | Trager | Log\[\(1\+x\^2\)/\(1\-x\^2\)\]/2 | Problem inverse\_function\_3 \(vi\) ArcTanh\(… | \{"failed", "transcendental function in the integrand"\} |
| product\_rule\-product\_rule\_2\-iv | Trager | \(a\*x\^5\-5\*a\*c\*x\^3\+2\*b\*x\^2\+2\*b\*c\)/\(\(x\^2\-c\)\^2\*Sqrt\[a\*x\^3\-b\]\) | Problem product\_rule\_2 \(iv\) \(solved, 2x … | \{"failed", "radicand is not a rational function over Q"\} |
| product\_rule\-product\_rule\_5\-v | Trager | \(9\*x\^9\-14\*b\*x\^6\-a\*x\^3\-4\*a\*b\)/\(x\^3\*\(x\^6\+a\)\^2\)\*\(\(x\^6\+a\)/\(x\^3\-b\)\)\^\(1/6\) | Problem product\_rule\_5 \(v\) | \{"failed", "radicand is not a rational function over Q"\} |
| realis\-prob\_realis\_1\-i | PMT | \(x\^2\-1\)/\(\(x\^2\-2\)\*\(x\^2\+1\)\*Sqrt\[3\-x\^4\]\) | Problem prob\_realis\_1 \(i\) | \{"needs torsion realisation \(milestone iii\)", \{\{\-2 \+ x\^2, \{\(1/6\*I\)/Sqrt\[2\], \(\-1/6\*I\)/Sqrt\[2\], \(\-1/6\*I\)/Sqrt\[2\]… |
| realis\-prob\_realis\_1\-iii | Trager | \(3\*b\*x\^4\+a\)/\(\(b\*c\*x\^6\+a\*c\*x\^2\+d\)\*Sqrt\[b\*x\^4\+a\]\) | Problem prob\_realis\_1 \(iii\) | \{"failed", "radicand is not a rational function over Q"\} |
| realis\-prob\_realis\_2\-iii\_upper | PMT | \(5\*x\-1\)/\(\(x\^2\-1\)\*Sqrt\[\(x\-1\)\^2\*\(x\+1\)\^3\-1\]\) | Problem prob\_realis\_2 \(iii\) upper signs | \{"needs torsion realisation \(milestone iii\)", \{\{\-1 \+ Sqrt\[x\], \{\-2\*I, 2\*I\}\}, \{1 \+ Sqrt\[x\], \{2\*I, \-2\*I\}\}, \{1 \+ x… |
| realis\-prob\_realis\_2\-iii\_lower | PMT | \(5\*x\+1\)/\(\(x\^2\-1\)\*Sqrt\[\(x\+1\)\^2\*\(x\-1\)\^3\-1\]\) | Problem prob\_realis\_2 \(iii\) lower signs | \{"needs torsion realisation \(milestone iii\)", \{\{\-1 \+ x, \{\-3\*I, 3\*I\}\}, \{1 \+ x, \{\-2\*I, 2\*I\}\}\}\} |
| realis\-prob\_realis\_2\-v\_upper\_n2 | PMT | \(5\*x\-1\)/\(\(x\^2\-1\)\*Sqrt\[\(x\+1\)\*\(x\^2\-1\)\^2\-1\]\) | Problem prob\_realis\_2 \(v\) n=2, numerator… | \{"needs torsion realisation \(milestone iii\)", \{\{\-1 \+ Sqrt\[x\], \{\-2\*I, 2\*I\}\}, \{1 \+ Sqrt\[x\], \{2\*I, \-2\*I\}\}, \{1 \+ x… |
| realis\-prob\_realis\_2\-v\_lower\_n2 | PMT | \(5\*x\+1\)/\(\(x\^2\-1\)\*Sqrt\[\(x\-1\)\*\(x\^2\-1\)\^2\-1\]\) | Problem prob\_realis\_2 \(v\) n=2, lower sig… | \{"needs torsion realisation \(milestone iii\)", \{\{\-1 \+ x, \{\-3\*I, 3\*I\}\}, \{1 \+ x, \{\-2\*I, 2\*I\}\}\}\} |
| realis\-prob\_realis\_2\-v\_upper\_n3 | PMT | \(7\*x\-1\)/\(\(x\^2\-1\)\*Sqrt\[\(x\+1\)\*\(x\^2\-1\)\^3\-1\]\) | Problem prob\_realis\_2 \(v\) n=3, numerator… | \{"needs torsion realisation \(milestone iii\)", \{\{\-1 \+ x, \{\-3\*I, 3\*I\}\}, \{1 \+ x, \{\-4\*I, 4\*I\}\}\}\} |
| realis\-prob\_realis\_2\-v\_lower\_n3 | PMT | \(7\*x\+1\)/\(\(x\^2\-1\)\*Sqrt\[\(x\-1\)\*\(x\^2\-1\)\^3\-1\]\) | Problem prob\_realis\_2 \(v\) n=3, lower sig… | \{"needs torsion realisation \(milestone iii\)", \{\{\-1 \+ Sqrt\[x\], \{\-4\*I, 4\*I\}\}, \{1 \+ Sqrt\[x\], \{4\*I, \-4\*I\}\}, \{1 \+ x… |
| realis\-prob\_realis\_3\-ii | Trager | \(3\*a\*x\^4\+2\*b\*x\^2\+c\)/\(\(a\*e\*x\^6\+b\*e\*x\^4\+c\*e\*x\^2\+d\)\*Sqrt\[a\*x\^4\+b\*x\^2\+c\]\) | Problem prob\_realis\_3 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| realis\-prob\_realis\_4\-i | PMT | 1/\(\(1\-x\^2\)\*Sqrt\[\(2\-x\^2\)\*\(x\^4\-2\*x\^2\+2\)\]\) | Problem prob\_realis\_4 \(i\) \(has a concise… | \{"needs torsion realisation \(milestone iii\)", \{\{\-1 \+ x, \{\-1/2, 1/2\}\}, \{1 \+ x, \{1/2, \-1/2\}\}\}\} |
| realis\-prob\_realis\_4\-ii | Trager | \(a\-2\*x\^2\)/\(\(2\*x\^4\-2\*a\*x\^2\+a\)\*Sqrt\[\(a\-x\^2\)\*\(x\^4\-a\*x\^2\+a\)\]\) | Problem prob\_realis\_4 \(ii\) \(a=2 gives \(i… | \{"failed", "radicand is not a rational function over Q"\} |
| nested\-intro\-1 | Trager | Sqrt\[a\*x\+b\*Sqrt\[x\]\] | Section intro, nested radicals \(closed f… | \{"failed", "nested radicals"\} |
| nested\-intro\-1\-inst | Trager | Sqrt\[2\*x\+3\*Sqrt\[x\]\] | Section intro, nested radicals \(closed f… | \{"failed", "nested radicals"\} |
| nested\-nested\_1\-i | Trager | Sqrt\[a/x\^2\+b\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(i\), n=\-2 | \{"failed", "nested radicals"\} |
| nested\-nested\_1\-i | PMT | Sqrt\[a/x\^2\+b\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(i\), n=\-2 | \{"needs torsion realisation \(milestone iii\)", \{\{Sqrt\[x\], \{2 Sqrt\[a\], \-2 Sqrt\[a\]\}\}\}\} |
| nested\-nested\_1\-i\-inst | Trager | Sqrt\[2/x\^2\+3\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(i\), n=\-2; inst a=2, b=… | \{"failed", "nested radicals"\} |
| nested\-nested\_1\-i\-inst | PMT | Sqrt\[2/x\^2\+3\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(i\), n=\-2; inst a=2, b=… | \{"needs torsion realisation \(milestone iii\)", \{\{Sqrt\[x\], \{2 Sqrt\[2\], \-2 Sqrt\[2\]\}\}\}\} |
| nested\-nested\_1\-ii | Trager | Sqrt\[a\*x\+b\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(ii\), n=1 | \{"failed", "nested radicals"\} |
| nested\-nested\_1\-ii | PMT | Sqrt\[a\*x\+b\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(ii\), n=1 | \{"failed", "verification failed"\} |
| nested\-nested\_1\-ii\-inst | Trager | Sqrt\[2\*x\+3\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(ii\), n=1; inst a=2, b=… | \{"failed", "nested radicals"\} |
| nested\-nested\_1\-iii | Trager | Sqrt\[a\*x\^2\+b\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(iii\), n=2 | \{"failed", "nested radicals"\} |
| nested\-nested\_1\-iii\-inst | Trager | Sqrt\[2\*x\^2\+3\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(iii\), n=2; inst a=2, b… | \{"failed", "nested radicals"\} |
| nested\-nested\_1\-iv | Trager | Sqrt\[a\*x\^3\+b\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(iv\), n=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_1\-iv\-inst | Trager | Sqrt\[2\*x\^3\+3\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(iv\), n=3; inst a=2, b=… | \{"failed", "nested radicals"\} |
| nested\-nested\_1\-v | Trager | Sqrt\[a\*x\^4\+b\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(v\), n=4 | \{"failed", "nested radicals"\} |
| nested\-nested\_1\-v | PMT | Sqrt\[a\*x\^4\+b\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(v\), n=4 | \{"failed", "no solution within bounds", \{12\}\} |
| nested\-nested\_1\-v\-inst | Trager | Sqrt\[2\*x\^4\+3\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(v\), n=4; inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_1\-v\-inst | PMT | Sqrt\[2\*x\^4\+3\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(v\), n=4; inst a=2, b=3 | \{"failed", "no solution within bounds", \{12\}\} |
| nested\-nested\_1\-vi | Trager | Sqrt\[a\*x\^7\+b\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(vi\), n=7 | \{"failed", "nested radicals"\} |
| nested\-nested\_1\-vi | PMT | Sqrt\[a\*x\^7\+b\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(vi\), n=7 | \{"failed", "no solution within bounds", \{18\}\} |
| nested\-nested\_1\-vi\-inst | Trager | Sqrt\[2\*x\^7\+3\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(vi\), n=7; inst a=2, b=… | \{"failed", "nested radicals"\} |
| nested\-nested\_1\-vi\-inst | PMT | Sqrt\[2\*x\^7\+3\*x\^2\*Sqrt\[x\]\] | Problem nested\_1 \(vi\), n=7; inst a=2, b=… | \{"failed", "no solution within bounds", \{18\}\} |
| nested\-nested\_2\-n1 | Trager | Sqrt\[a\*x\+b\*x\^2/Sqrt\[x\]\] | Problem nested\_2, n=1 \(algebraic for all… | \{"failed", "nested radicals"\} |
| nested\-nested\_2\-n1\-inst | Trager | Sqrt\[2\*x\+3\*x\^2/Sqrt\[x\]\] | Problem nested\_2, n=1 \(algebraic for all… | \{"failed", "nested radicals"\} |
| nested\-nested\_2\-n2 | Trager | Sqrt\[a\*x\^2\+b\*x\^3/Sqrt\[x\]\] | Problem nested\_2, n=2 \(algebraic for all… | \{"failed", "nested radicals"\} |
| nested\-nested\_2\-n2\-inst | Trager | Sqrt\[2\*x\^2\+3\*x\^3/Sqrt\[x\]\] | Problem nested\_2, n=2 \(algebraic for all… | \{"failed", "nested radicals"\} |
| nested\-nested\_2\-n3 | Trager | Sqrt\[a\*x\^3\+b\*x\^4/Sqrt\[x\]\] | Problem nested\_2, n=3 \(algebraic for all… | \{"failed", "nested radicals"\} |
| nested\-nested\_2\-n3\-inst | Trager | Sqrt\[2\*x\^3\+3\*x\^4/Sqrt\[x\]\] | Problem nested\_2, n=3 \(algebraic for all… | \{"failed", "nested radicals"\} |
| nested\-nested\_3\-i | Trager | \(a/x\^3\+b\*x\*Sqrt\[x\]\)\^\(1/3\) | Problem nested\_3 \(i\), n=\-3 | \{"failed", "nested radicals"\} |
| nested\-nested\_3\-i\-inst | Trager | \(2/x\^3\+3\*x\*Sqrt\[x\]\)\^\(1/3\) | Problem nested\_3 \(i\), n=\-3; inst a=2, b=… | \{"failed", "nested radicals"\} |
| nested\-nested\_3\-ii | Trager | \(a\+b\*x\*Sqrt\[x\]\)\^\(1/3\) | Problem nested\_3 \(ii\), n=0 | \{"failed", "nested radicals"\} |
| nested\-nested\_3\-ii | PMT | \(a\+b\*x\*Sqrt\[x\]\)\^\(1/3\) | Problem nested\_3 \(ii\), n=0 | \{"failed", "internal: solution does not verify"\} |
| nested\-nested\_3\-ii\-inst | Trager | \(2\+3\*x\*Sqrt\[x\]\)\^\(1/3\) | Problem nested\_3 \(ii\), n=0; inst a=2, b=… | \{"failed", "nested radicals"\} |
| nested\-nested\_3\-iii | Trager | \(a\*x\+b\*x\*Sqrt\[x\]\)\^\(1/3\) | Problem nested\_3 \(iii\), n=1 | \{"failed", "nested radicals"\} |
| nested\-nested\_3\-iii\-inst | Trager | \(2\*x\+3\*x\*Sqrt\[x\]\)\^\(1/3\) | Problem nested\_3 \(iii\), n=1; inst a=2, b… | \{"failed", "nested radicals"\} |
| nested\-nested\_3\-iv | Trager | \(a\*x\^2\+b\*x\*Sqrt\[x\]\)\^\(1/3\) | Problem nested\_3 \(iv\), n=2 | \{"failed", "nested radicals"\} |
| nested\-nested\_3\-iv\-inst | Trager | \(2\*x\^2\+3\*x\*Sqrt\[x\]\)\^\(1/3\) | Problem nested\_3 \(iv\), n=2; inst a=2, b=… | \{"failed", "nested radicals"\} |
| nested\-nested\_3\-v | Trager | \(a\*x\^3\+b\*x\*Sqrt\[x\]\)\^\(1/3\) | Problem nested\_3 \(v\), n=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_3\-v | PMT | \(a\*x\^3\+b\*x\*Sqrt\[x\]\)\^\(1/3\) | Problem nested\_3 \(v\), n=3 | \{"failed", "internal: solution does not verify"\} |
| nested\-nested\_3\-v\-inst | Trager | \(2\*x\^3\+3\*x\*Sqrt\[x\]\)\^\(1/3\) | Problem nested\_3 \(v\), n=3; inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_4\-ex1 | Trager | \(a/x\^4\+b\*x\^\(4/3\)\)\^\(1/4\) | Problem nested\_4, solution example 1 \(n=… | \{"failed", "nested radicals"\} |
| nested\-nested\_4\-ex1\-inst | Trager | \(2/x\^4\+3\*x\^\(4/3\)\)\^\(1/4\) | Problem nested\_4, solution example 1 \(n=… | \{"failed", "nested radicals"\} |
| nested\-nested\_4\-ex2 | Trager | \(a\*x\^3\+b\*x\^\(16/3\)\)\^\(1/4\) | Problem nested\_4, solution example 2 \(n=… | \{"failed", "nested radicals"\} |
| nested\-nested\_4\-ex2 | PMT | \(a\*x\^3\+b\*x\^\(16/3\)\)\^\(1/4\) | Problem nested\_4, solution example 2 \(n=… | \{"failed", "no solution within bounds", \{28\}\} |
| nested\-nested\_4\-ex2\-inst | Trager | \(2\*x\^3\+3\*x\^\(16/3\)\)\^\(1/4\) | Problem nested\_4, solution example 2 \(n=… | \{"failed", "nested radicals"\} |
| nested\-nested\_4\-ex2\-inst | PMT | \(2\*x\^3\+3\*x\^\(16/3\)\)\^\(1/4\) | Problem nested\_4, solution example 2 \(n=… | \{"failed", "no solution within bounds", \{28\}\} |
| nested\-nested\_5\-i | Trager | Sqrt\[d\*x\+c\*x\^2\*Sqrt\[a/x\+b\]\]/Sqrt\[a/x\+b\] | Problem nested\_5 \(i\) | \{"failed", "nested radicals"\} |
| nested\-nested\_5\-i | PMT | Sqrt\[d\*x\+c\*x\^2\*Sqrt\[a/x\+b\]\]/Sqrt\[a/x\+b\] | Problem nested\_5 \(i\) | \{"failed", "tower construction failed", Sqrt\[d x \+ c Sqrt\[b \+ a/x\] x\^2\]/Sqrt\[b \+ a/x\]\} |
| nested\-nested\_5\-i\-inst | Trager | Sqrt\[7\*x\+5\*x\^2\*Sqrt\[2/x\+3\]\]/Sqrt\[2/x\+3\] | Problem nested\_5 \(i\); inst a=2, b=3, c=5… | \{"failed", "nested radicals"\} |
| nested\-nested\_5\-i\-inst | PMT | Sqrt\[7\*x\+5\*x\^2\*Sqrt\[2/x\+3\]\]/Sqrt\[2/x\+3\] | Problem nested\_5 \(i\); inst a=2, b=3, c=5… | \{"failed", "tower construction failed", Sqrt\[7 x \+ 5 Sqrt\[3 \+ 2/x\] x\^2\]/Sqrt\[3 \+ 2/x\]\} |
| nested\-nested\_5\-ii | Trager | Sqrt\[a/x\+b\]/Sqrt\[d\*x\+c\*x\^2\*Sqrt\[a/x\+b\]\] | Problem nested\_5 \(ii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_5\-ii | PMT | Sqrt\[a/x\+b\]/Sqrt\[d\*x\+c\*x\^2\*Sqrt\[a/x\+b\]\] | Problem nested\_5 \(ii\) | \{"failed", "tower construction failed", Sqrt\[b \+ a/x\]/Sqrt\[d x \+ c Sqrt\[b \+ a/x\] x\^2\]\} |
| nested\-nested\_5\-ii\-inst | Trager | Sqrt\[2/x\+3\]/Sqrt\[7\*x\+5\*x\^2\*Sqrt\[2/x\+3\]\] | Problem nested\_5 \(ii\); inst a=2, b=3, c=… | \{"failed", "nested radicals"\} |
| nested\-nested\_5\-ii\-inst | PMT | Sqrt\[2/x\+3\]/Sqrt\[7\*x\+5\*x\^2\*Sqrt\[2/x\+3\]\] | Problem nested\_5 \(ii\); inst a=2, b=3, c=… | \{"failed", "tower construction failed", Sqrt\[3 \+ 2/x\]/Sqrt\[7 x \+ 5 Sqrt\[3 \+ 2/x\] x\^2\]\} |
| nested\-nested\_5\-iii | Trager | Sqrt\[d\*x\+c\*x\^2\*Sqrt\[a/x\+b\]\]/\(\(c\*x\^2\+d\*x\)\*Sqrt\[a/x\+b\]\) | Problem nested\_5 \(iii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_5\-iii | PMT | Sqrt\[d\*x\+c\*x\^2\*Sqrt\[a/x\+b\]\]/\(\(c\*x\^2\+d\*x\)\*Sqrt\[a/x\+b\]\) | Problem nested\_5 \(iii\) | \{"failed", "tower construction failed", Sqrt\[d x \+ c Sqrt\[b \+ a/x\] x\^2\]/\(Sqrt\[b \+ a/x\] \(d x \+ c x\^2\)\)\} |
| nested\-nested\_5\-iii\-inst | Trager | Sqrt\[7\*x\+5\*x\^2\*Sqrt\[2/x\+3\]\]/\(\(5\*x\^2\+7\*x\)\*Sqrt\[2/x\+3\]\) | Problem nested\_5 \(iii\); inst a=2, b=3, c… | \{"failed", "nested radicals"\} |
| nested\-nested\_5\-iii\-inst | PMT | Sqrt\[7\*x\+5\*x\^2\*Sqrt\[2/x\+3\]\]/\(\(5\*x\^2\+7\*x\)\*Sqrt\[2/x\+3\]\) | Problem nested\_5 \(iii\); inst a=2, b=3, c… | \{"failed", "tower construction failed", Sqrt\[7 x \+ 5 Sqrt\[3 \+ 2/x\] x\^2\]/\(Sqrt\[3 \+ 2/x\] \(7 x \+ 5 x\^2\)\)\} |
| nested\-nested\_5\-iv | Trager | Sqrt\[a/x\+b\]/\(\(c\*x\^2\+d\*x\)\*Sqrt\[d\*x\+c\*x\^2\*Sqrt\[a/x\+b\]\]\) | Problem nested\_5 \(iv\) | \{"failed", "nested radicals"\} |
| nested\-nested\_5\-iv | PMT | Sqrt\[a/x\+b\]/\(\(c\*x\^2\+d\*x\)\*Sqrt\[d\*x\+c\*x\^2\*Sqrt\[a/x\+b\]\]\) | Problem nested\_5 \(iv\) | \{"failed", "tower construction failed", Sqrt\[b \+ a/x\]/\(\(d x \+ c x\^2\) Sqrt\[d x \+ c Sqrt\[b \+ a/x\] x\^2\]\)\} |
| nested\-nested\_5\-iv\-inst | Trager | Sqrt\[2/x\+3\]/\(\(5\*x\^2\+7\*x\)\*Sqrt\[7\*x\+5\*x\^2\*Sqrt\[2/x\+3\]\]\) | Problem nested\_5 \(iv\); inst a=2, b=3, c=… | \{"failed", "nested radicals"\} |
| nested\-nested\_5\-iv\-inst | PMT | Sqrt\[2/x\+3\]/\(\(5\*x\^2\+7\*x\)\*Sqrt\[7\*x\+5\*x\^2\*Sqrt\[2/x\+3\]\]\) | Problem nested\_5 \(iv\); inst a=2, b=3, c=… | \{"failed", "tower construction failed", Sqrt\[3 \+ 2/x\]/\(\(7 x \+ 5 x\^2\) Sqrt\[7 x \+ 5 Sqrt\[3 \+ 2/x\] x\^2\]\)\} |
| nested\-nested\_6\-i | Trager | Sqrt\[d\*x\^2\+c\*x\^2/Sqrt\[a/x\+b\]\]/Sqrt\[a/x\+b\] | Problem nested\_6 \(i\) | \{"failed", "nested radicals"\} |
| nested\-nested\_6\-i\-inst | Trager | Sqrt\[7\*x\^2\+5\*x\^2/Sqrt\[2/x\+3\]\]/Sqrt\[2/x\+3\] | Problem nested\_6 \(i\); inst a=2, b=3, c=5… | \{"failed", "nested radicals"\} |
| nested\-nested\_6\-ii | Trager | Sqrt\[a/x\+b\]/Sqrt\[d\*x\^2\+c\*x\^2/Sqrt\[a/x\+b\]\] | Problem nested\_6 \(ii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_6\-ii\-inst | Trager | Sqrt\[2/x\+3\]/Sqrt\[7\*x\^2\+5\*x\^2/Sqrt\[2/x\+3\]\] | Problem nested\_6 \(ii\); inst a=2, b=3, c=… | \{"failed", "nested radicals"\} |
| nested\-nested\_6\-iii | Trager | Sqrt\[a/x\+b\]\*Sqrt\[d\*x\^2\+c\*x\^2/Sqrt\[a/x\+b\]\] | Problem nested\_6 \(iii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_6\-iii\-inst | Trager | Sqrt\[2/x\+3\]\*Sqrt\[7\*x\^2\+5\*x\^2/Sqrt\[2/x\+3\]\] | Problem nested\_6 \(iii\); inst a=2, b=3, c… | \{"failed", "nested radicals"\} |
| nested\-nested\_6\-iv | Trager | 1/\(Sqrt\[a/x\+b\]\*Sqrt\[d\*x\^2\+c\*x\^2/Sqrt\[a/x\+b\]\]\) | Problem nested\_6 \(iv\) | \{"failed", "nested radicals"\} |
| nested\-nested\_6\-iv\-inst | Trager | 1/\(Sqrt\[2/x\+3\]\*Sqrt\[7\*x\^2\+5\*x\^2/Sqrt\[2/x\+3\]\]\) | Problem nested\_6 \(iv\); inst a=2, b=3, c=… | \{"failed", "nested radicals"\} |
| nested\-nested\_7\-i | Trager | Sqrt\[a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\]/\(\(b\*x\+c\)\*Sqrt\[a\^2\*x\^2\+b\*x\+c\]\) | Problem nested\_7 \(i\) | \{"failed", "nested radicals"\} |
| nested\-nested\_7\-i | PMT | Sqrt\[a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\]/\(\(b\*x\+c\)\*Sqrt\[a\^2\*x\^2\+b\*x\+c\]\) | Problem nested\_7 \(i\) | \{"failed", "tower construction failed", Sqrt\[a x \+ Sqrt\[c \+ b x \+ a\^2 x\^2\]\]/\(\(c \+ b x\) Sqrt\[c \+ b x \+ a\^2 x\^2\]… |
| nested\-nested\_7\-i\-inst | Trager | Sqrt\[2\*x\+Sqrt\[2\^2\*x\^2\+3\*x\+5\]\]/\(\(3\*x\+5\)\*Sqrt\[2\^2\*x\^2\+3\*x\+5\]\) | Problem nested\_7 \(i\); inst a=2, b=3, c=5 | \{"failed", "nested radicals"\} |
| nested\-nested\_7\-ii | Trager | \(a\^2\*x\^2\+b\*x\+c\)/Sqrt\[a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\] | Problem nested\_7 \(ii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_7\-ii | PMT | \(a\^2\*x\^2\+b\*x\+c\)/Sqrt\[a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\] | Problem nested\_7 \(ii\) | \{"failed", "tower construction failed", \(c \+ b x \+ a\^2 x\^2\)/Sqrt\[a x \+ Sqrt\[c \+ b x \+ a\^2 x\^2\]\]\} |
| nested\-nested\_7\-ii\-inst | Trager | \(2\^2\*x\^2\+3\*x\+5\)/Sqrt\[2\*x\+Sqrt\[2\^2\*x\^2\+3\*x\+5\]\] | Problem nested\_7 \(ii\); inst a=2, b=3, c=… | \{"failed", "nested radicals"\} |
| nested\-nested\_7\-iii | Trager | Sqrt\[a\^2\*x\^2\+b\*x\+c\]/Sqrt\[a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\] | Problem nested\_7 \(iii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_7\-iii | PMT | Sqrt\[a\^2\*x\^2\+b\*x\+c\]/Sqrt\[a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\] | Problem nested\_7 \(iii\) | \{"failed", "tower construction failed", Sqrt\[c \+ b x \+ a\^2 x\^2\]/Sqrt\[a x \+ Sqrt\[c \+ b x \+ a\^2 x\^2\]\]\} |
| nested\-nested\_7\-iii\-inst | Trager | Sqrt\[2\^2\*x\^2\+3\*x\+5\]/Sqrt\[2\*x\+Sqrt\[2\^2\*x\^2\+3\*x\+5\]\] | Problem nested\_7 \(iii\); inst a=2, b=3, c… | \{"failed", "nested radicals"\} |
| nested\-nested\_7\-iv | Trager | \(a\*x\^2\-b\*x\-c\)/Sqrt\[a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\] | Problem nested\_7 \(iv\) | \{"failed", "nested radicals"\} |
| nested\-nested\_7\-iv | PMT | \(a\*x\^2\-b\*x\-c\)/Sqrt\[a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\] | Problem nested\_7 \(iv\) | \{"failed", "tower construction failed", \(\-c \- b x \+ a x\^2\)/Sqrt\[a x \+ Sqrt\[c \+ b x \+ a\^2 x\^2\]\]\} |
| nested\-nested\_7\-iv\-inst | Trager | \(2\*x\^2\-3\*x\-5\)/Sqrt\[2\*x\+Sqrt\[2\^2\*x\^2\+3\*x\+5\]\] | Problem nested\_7 \(iv\); inst a=2, b=3, c=… | \{"failed", "nested radicals"\} |
| nested\-nested\_7\-v | Trager | Sqrt\[a\^2\*x\^2\+b\*x\+c\]/\(x\^2\*Sqrt\[\-a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\]\) | Problem nested\_7 \(v\) | \{"failed", "nested radicals"\} |
| nested\-nested\_7\-v | PMT | Sqrt\[a\^2\*x\^2\+b\*x\+c\]/\(x\^2\*Sqrt\[\-a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\]\) | Problem nested\_7 \(v\) | \{"failed", "tower construction failed", Sqrt\[c \+ b x \+ a\^2 x\^2\]/\(x\^2 Sqrt\[\-a x \+ Sqrt\[c \+ b x \+ a\^2 x\^2\]\]\)\} |
| nested\-nested\_7\-v\-inst | Trager | Sqrt\[2\^2\*x\^2\+3\*x\+5\]/\(x\^2\*Sqrt\[\-2\*x\+Sqrt\[2\^2\*x\^2\+3\*x\+5\]\]\) | Problem nested\_7 \(v\); inst a=2, b=3, c=5 | \{"failed", "nested radicals"\} |
| nested\-nested\_8\-i | Trager | Sqrt\[d\+b\*Sqrt\[x\]\+Sqrt\[a\*Sqrt\[x\]\+b\^2\*x\+c\]\]/x | Problem nested\_8 \(i\) | \{"failed", "nested radicals"\} |
| nested\-nested\_8\-i | PMT | Sqrt\[d\+b\*Sqrt\[x\]\+Sqrt\[a\*Sqrt\[x\]\+b\^2\*x\+c\]\]/x | Problem nested\_8 \(i\) | \{"failed", "tower construction failed", Sqrt\[d \+ b Sqrt\[x\] \+ Sqrt\[c \+ a Sqrt\[x\] \+ b\^2 x\]\]/x\} |
| nested\-nested\_8\-i\-inst | Trager | Sqrt\[7\+3\*Sqrt\[x\]\+Sqrt\[2\*Sqrt\[x\]\+3\^2\*x\+5\]\]/x | Problem nested\_8 \(i\); inst a=2, b=3, c=5… | \{"failed", "nested radicals"\} |
| nested\-nested\_8\-ii | Trager | Sqrt\[d\+b\*Sqrt\[x\]\+Sqrt\[a\*Sqrt\[x\]\+b\^2\*x\+c\]\] | Problem nested\_8 \(ii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_8\-ii | PMT | Sqrt\[d\+b\*Sqrt\[x\]\+Sqrt\[a\*Sqrt\[x\]\+b\^2\*x\+c\]\] | Problem nested\_8 \(ii\) | \{"failed", "tower construction failed", Sqrt\[d \+ b Sqrt\[x\] \+ Sqrt\[c \+ a Sqrt\[x\] \+ b\^2 x\]\]\} |
| nested\-nested\_8\-ii\-inst | Trager | Sqrt\[7\+3\*Sqrt\[x\]\+Sqrt\[2\*Sqrt\[x\]\+3\^2\*x\+5\]\] | Problem nested\_8 \(ii\); inst a=2, b=3, c=… | \{"failed", "nested radicals"\} |
| nested\-nested\_9\-i | Trager | 1/\(\-a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\)\^\(1/4\) | Problem nested\_9 \(i\) | \{"failed", "nested radicals"\} |
| nested\-nested\_9\-i | PMT | 1/\(\-a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\)\^\(1/4\) | Problem nested\_9 \(i\) | \{"failed", "tower construction failed", 1/\(\-a x \+ Sqrt\[c \+ b x \+ a\^2 x\^2\]\)\^\(1/4\)\} |
| nested\-nested\_9\-i\-inst | Trager | 1/\(\-2\*x\+Sqrt\[2\^2\*x\^2\+3\*x\+5\]\)\^\(1/4\) | Problem nested\_9 \(i\); inst a=2, b=3, c=5 | \{"failed", "nested radicals"\} |
| nested\-nested\_9\-ii | Trager | 1/\(\-a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\)\^\(3/4\) | Problem nested\_9 \(ii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_9\-ii | PMT | 1/\(\-a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\)\^\(3/4\) | Problem nested\_9 \(ii\) | \{"failed", "tower construction failed", 1/\(\-a x \+ Sqrt\[c \+ b x \+ a\^2 x\^2\]\)\^\(3/4\)\} |
| nested\-nested\_9\-ii\-inst | Trager | 1/\(\-2\*x\+Sqrt\[2\^2\*x\^2\+3\*x\+5\]\)\^\(3/4\) | Problem nested\_9 \(ii\); inst a=2, b=3, c=… | \{"failed", "nested radicals"\} |
| nested\-nested\_9\-iii | Trager | \(\-a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\)\^\(1/4\)/\(b\*x\+c\) | Problem nested\_9 \(iii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_9\-iii | PMT | \(\-a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\)\^\(1/4\)/\(b\*x\+c\) | Problem nested\_9 \(iii\) | \{"failed", "tower construction failed", \(\-a x \+ Sqrt\[c \+ b x \+ a\^2 x\^2\]\)\^\(1/4\)/\(c \+ b x\)\} |
| nested\-nested\_9\-iii\-inst | Trager | \(\-2\*x\+Sqrt\[2\^2\*x\^2\+3\*x\+5\]\)\^\(1/4\)/\(3\*x\+5\) | Problem nested\_9 \(iii\); inst a=2, b=3, c… | \{"failed", "nested radicals"\} |
| nested\-nested\_9\-iv | Trager | 1/\(\-a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\)\^\(5/6\) | Problem nested\_9 \(iv\) | \{"failed", "nested radicals"\} |
| nested\-nested\_9\-iv | PMT | 1/\(\-a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\)\^\(5/6\) | Problem nested\_9 \(iv\) | \{"failed", "tower construction failed", 1/\(\-a x \+ Sqrt\[c \+ b x \+ a\^2 x\^2\]\)\^\(5/6\)\} |
| nested\-nested\_9\-iv\-inst | Trager | 1/\(\-2\*x\+Sqrt\[2\^2\*x\^2\+3\*x\+5\]\)\^\(5/6\) | Problem nested\_9 \(iv\); inst a=2, b=3, c=… | \{"failed", "nested radicals"\} |
| nested\-nested\_9\-v | Trager | x/\(\-a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\)\^\(1/8\) | Problem nested\_9 \(v\) | \{"failed", "nested radicals"\} |
| nested\-nested\_9\-v | PMT | x/\(\-a\*x\+Sqrt\[a\^2\*x\^2\+b\*x\+c\]\)\^\(1/8\) | Problem nested\_9 \(v\) | \{"failed", "tower construction failed", x/\(\-a x \+ Sqrt\[c \+ b x \+ a\^2 x\^2\]\)\^\(1/8\)\} |
| nested\-nested\_9\-v\-inst | Trager | x/\(\-2\*x\+Sqrt\[2\^2\*x\^2\+3\*x\+5\]\)\^\(1/8\) | Problem nested\_9 \(v\); inst a=2, b=3, c=5 | \{"failed", "nested radicals"\} |
| nested\-nested\_10\-i | Trager | Sqrt\[a\^2\*x\^2\+b\]/Sqrt\[a\*x\^2\+x\*Sqrt\[a\^2\*x\^2\+b\]\] | Problem nested\_10 \(i\) | \{"failed", "nested radicals"\} |
| nested\-nested\_10\-i | PMT | Sqrt\[a\^2\*x\^2\+b\]/Sqrt\[a\*x\^2\+x\*Sqrt\[a\^2\*x\^2\+b\]\] | Problem nested\_10 \(i\) | \{"failed", "tower construction failed", Sqrt\[b \+ a\^2 x\^2\]/Sqrt\[a x\^2 \+ x Sqrt\[b \+ a\^2 x\^2\]\]\} |
| nested\-nested\_10\-i\-inst | Trager | Sqrt\[2\^2\*x\^2\+3\]/Sqrt\[2\*x\^2\+x\*Sqrt\[2\^2\*x\^2\+3\]\] | Problem nested\_10 \(i\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_10\-i\-inst | PMT | Sqrt\[2\^2\*x\^2\+3\]/Sqrt\[2\*x\^2\+x\*Sqrt\[2\^2\*x\^2\+3\]\] | Problem nested\_10 \(i\); inst a=2, b=3 | \{"failed", "tower construction failed", Sqrt\[3 \+ 4 x\^2\]/Sqrt\[2 x\^2 \+ x Sqrt\[3 \+ 4 x\^2\]\]\} |
| nested\-nested\_10\-ii | Trager | Sqrt\[a\^2\*x\^2\+b\]\*Sqrt\[a\*x\^2\+x\*Sqrt\[a\^2\*x\^2\+b\]\]/x | Problem nested\_10 \(ii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_10\-ii | PMT | Sqrt\[a\^2\*x\^2\+b\]\*Sqrt\[a\*x\^2\+x\*Sqrt\[a\^2\*x\^2\+b\]\]/x | Problem nested\_10 \(ii\) | \{"failed", "tower construction failed", \(Sqrt\[b \+ a\^2 x\^2\] Sqrt\[a x\^2 \+ x Sqrt\[b \+ a\^2 x\^2\]\]\)/x\} |
| nested\-nested\_10\-ii\-inst | Trager | Sqrt\[2\^2\*x\^2\+3\]\*Sqrt\[2\*x\^2\+x\*Sqrt\[2\^2\*x\^2\+3\]\]/x | Problem nested\_10 \(ii\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_10\-ii\-inst | PMT | Sqrt\[2\^2\*x\^2\+3\]\*Sqrt\[2\*x\^2\+x\*Sqrt\[2\^2\*x\^2\+3\]\]/x | Problem nested\_10 \(ii\); inst a=2, b=3 | \{"failed", "tower construction failed", \(Sqrt\[3 \+ 4 x\^2\] Sqrt\[2 x\^2 \+ x Sqrt\[3 \+ 4 x\^2\]\]\)/x\} |
| nested\-nested\_10\-iii | Trager | x\^2/\(Sqrt\[a\^2\*x\^2\+b\]\*Sqrt\[\-a\*x\^2\+x\*Sqrt\[a\^2\*x\^2\+b\]\]\) | Problem nested\_10 \(iii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_10\-iii | PMT | x\^2/\(Sqrt\[a\^2\*x\^2\+b\]\*Sqrt\[\-a\*x\^2\+x\*Sqrt\[a\^2\*x\^2\+b\]\]\) | Problem nested\_10 \(iii\) | \{"failed", "tower construction failed", x\^2/\(Sqrt\[b \+ a\^2 x\^2\] Sqrt\[\-a x\^2 \+ x Sqrt\[b \+ a\^2 x\^2\]\]\)\} |
| nested\-nested\_10\-iii\-inst | Trager | x\^2/\(Sqrt\[2\^2\*x\^2\+3\]\*Sqrt\[\-2\*x\^2\+x\*Sqrt\[2\^2\*x\^2\+3\]\]\) | Problem nested\_10 \(iii\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_10\-iii\-inst | PMT | x\^2/\(Sqrt\[2\^2\*x\^2\+3\]\*Sqrt\[\-2\*x\^2\+x\*Sqrt\[2\^2\*x\^2\+3\]\]\) | Problem nested\_10 \(iii\); inst a=2, b=3 | \{"failed", "tower construction failed", x\^2/\(Sqrt\[3 \+ 4 x\^2\] Sqrt\[\-2 x\^2 \+ x Sqrt\[3 \+ 4 x\^2\]\]\)\} |
| nested\-nested\_10\-iv | Trager | Sqrt\[\-a\*x\^2\+x\*Sqrt\[a\^2\*x\^2\+b\]\]/\(a\^2\*x\^2\+b\*x\) | Problem nested\_10 \(iv\) | \{"failed", "nested radicals"\} |
| nested\-nested\_10\-iv | PMT | Sqrt\[\-a\*x\^2\+x\*Sqrt\[a\^2\*x\^2\+b\]\]/\(a\^2\*x\^2\+b\*x\) | Problem nested\_10 \(iv\) | \{"failed", "tower construction failed", Sqrt\[\-a x\^2 \+ x Sqrt\[b \+ a\^2 x\^2\]\]/\(b x \+ a\^2 x\^2\)\} |
| nested\-nested\_10\-iv\-inst | Trager | Sqrt\[\-2\*x\^2\+x\*Sqrt\[2\^2\*x\^2\+3\]\]/\(2\^2\*x\^2\+3\*x\) | Problem nested\_10 \(iv\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_10\-iv\-inst | PMT | Sqrt\[\-2\*x\^2\+x\*Sqrt\[2\^2\*x\^2\+3\]\]/\(2\^2\*x\^2\+3\*x\) | Problem nested\_10 \(iv\); inst a=2, b=3 | \{"failed", "tower construction failed", Sqrt\[\-2 x\^2 \+ x Sqrt\[3 \+ 4 x\^2\]\]/\(3 x \+ 4 x\^2\)\} |
| nested\-nested\_10\-v | Trager | \(a\^2\*x\^2\-b\)\^2/\(\(a\^2\*x\^2\+b\)\^2\*Sqrt\[a\*x\^2\+x\*Sqrt\[a\^2\*x\^2\+b\]\]\) | Problem nested\_10 \(v\) | \{"failed", "nested radicals"\} |
| nested\-nested\_10\-v | PMT | \(a\^2\*x\^2\-b\)\^2/\(\(a\^2\*x\^2\+b\)\^2\*Sqrt\[a\*x\^2\+x\*Sqrt\[a\^2\*x\^2\+b\]\]\) | Problem nested\_10 \(v\) | \{"failed", "tower construction failed", \(\-b \+ a\^2 x\^2\)\^2/\(\(b \+ a\^2 x\^2\)\^2 Sqrt\[a x\^2 \+ x Sqrt\[b \+ a\^2 x\^2\]\]\)\} |
| nested\-nested\_10\-v\-inst | Trager | \(2\^2\*x\^2\-3\)\^2/\(\(2\^2\*x\^2\+3\)\^2\*Sqrt\[2\*x\^2\+x\*Sqrt\[2\^2\*x\^2\+3\]\]\) | Problem nested\_10 \(v\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_10\-v\-inst | PMT | \(2\^2\*x\^2\-3\)\^2/\(\(2\^2\*x\^2\+3\)\^2\*Sqrt\[2\*x\^2\+x\*Sqrt\[2\^2\*x\^2\+3\]\]\) | Problem nested\_10 \(v\); inst a=2, b=3 | \{"failed", "tower construction failed", \(\-3 \+ 4 x\^2\)\^2/\(\(3 \+ 4 x\^2\)\^2 Sqrt\[2 x\^2 \+ x Sqrt\[3 \+ 4 x\^2\]\]\)\} |
| nested\-nested\_11\-i | Trager | Sqrt\[x\^2\+Sqrt\[x\^4\+1\]\]/\(\(x\+1\)\^2\*Sqrt\[x\^4\+1\]\) | Problem nested\_11 \(i\) | \{"failed", "nested radicals"\} |
| nested\-nested\_11\-i | PMT | Sqrt\[x\^2\+Sqrt\[x\^4\+1\]\]/\(\(x\+1\)\^2\*Sqrt\[x\^4\+1\]\) | Problem nested\_11 \(i\) | \{"failed", "tower construction failed", Sqrt\[x\^2 \+ Sqrt\[1 \+ x\^4\]\]/\(\(1 \+ x\)\^2 Sqrt\[1 \+ x\^4\]\)\} |
| nested\-nested\_11\-ii | Trager | Sqrt\[a\*x\^2\+Sqrt\[a\^2\*x\^4\+b\]\]/\(\(c\*x\+d\)\^2\*Sqrt\[a\^2\*x\^4\+b\]\) | Problem nested\_11 \(ii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_11\-ii | PMT | Sqrt\[a\*x\^2\+Sqrt\[a\^2\*x\^4\+b\]\]/\(\(c\*x\+d\)\^2\*Sqrt\[a\^2\*x\^4\+b\]\) | Problem nested\_11 \(ii\) | \{"failed", "tower construction failed", Sqrt\[a x\^2 \+ Sqrt\[b \+ a\^2 x\^4\]\]/\(\(d \+ c x\)\^2 Sqrt\[b \+ a\^2 x\^4\]\)\} |
| nested\-nested\_11\-ii\-inst | Trager | Sqrt\[2\*x\^2\+Sqrt\[2\^2\*x\^4\+3\]\]/\(\(5\*x\+7\)\^2\*Sqrt\[2\^2\*x\^4\+3\]\) | Problem nested\_11 \(ii\); inst a=2, b=3, c… | \{"failed", "nested radicals"\} |
| nested\-nested\_11\-ii\-inst | PMT | Sqrt\[2\*x\^2\+Sqrt\[2\^2\*x\^4\+3\]\]/\(\(5\*x\+7\)\^2\*Sqrt\[2\^2\*x\^4\+3\]\) | Problem nested\_11 \(ii\); inst a=2, b=3, c… | \{"failed", "tower construction failed", Sqrt\[2 x\^2 \+ Sqrt\[3 \+ 4 x\^4\]\]/\(\(7 \+ 5 x\)\^2 Sqrt\[3 \+ 4 x\^4\]\)\} |
| nested\-nested\_11\-iii | Trager | Sqrt\[a\*x\^2\+Sqrt\[a\^2\*x\^4\+b\]\] | Problem nested\_11 \(iii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_11\-iii | PMT | Sqrt\[a\*x\^2\+Sqrt\[a\^2\*x\^4\+b\]\] | Problem nested\_11 \(iii\) | \{"failed", "tower construction failed", Sqrt\[a x\^2 \+ Sqrt\[b \+ a\^2 x\^4\]\]\} |
| nested\-nested\_11\-iii\-inst | Trager | Sqrt\[2\*x\^2\+Sqrt\[2\^2\*x\^4\+3\]\] | Problem nested\_11 \(iii\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_11\-iii\-inst | PMT | Sqrt\[2\*x\^2\+Sqrt\[2\^2\*x\^4\+3\]\] | Problem nested\_11 \(iii\); inst a=2, b=3 | \{"failed", "tower construction failed", Sqrt\[2 x\^2 \+ Sqrt\[3 \+ 4 x\^4\]\]\} |
| nested\-nested\_11\-iv | Trager | \(a\^2\*x\^4\+b\)\^2/\(\(a\^2\*x\^4\-b\)\^2\*Sqrt\[a\*x\^2\+Sqrt\[a\^2\*x\^4\+b\]\]\) | Problem nested\_11 \(iv\) | \{"failed", "nested radicals"\} |
| nested\-nested\_11\-iv | PMT | \(a\^2\*x\^4\+b\)\^2/\(\(a\^2\*x\^4\-b\)\^2\*Sqrt\[a\*x\^2\+Sqrt\[a\^2\*x\^4\+b\]\]\) | Problem nested\_11 \(iv\) | \{"failed", "tower construction failed", \(b \+ a\^2 x\^4\)\^2/\(\(\-b \+ a\^2 x\^4\)\^2 Sqrt\[a x\^2 \+ Sqrt\[b \+ a\^2 x\^4\]\]\)\} |
| nested\-nested\_11\-iv\-inst | Trager | \(2\^2\*x\^4\+3\)\^2/\(\(2\^2\*x\^4\-3\)\^2\*Sqrt\[2\*x\^2\+Sqrt\[2\^2\*x\^4\+3\]\]\) | Problem nested\_11 \(iv\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_11\-iv\-inst | PMT | \(2\^2\*x\^4\+3\)\^2/\(\(2\^2\*x\^4\-3\)\^2\*Sqrt\[2\*x\^2\+Sqrt\[2\^2\*x\^4\+3\]\]\) | Problem nested\_11 \(iv\); inst a=2, b=3 | \{"failed", "tower construction failed", \(3 \+ 4 x\^4\)\^2/\(\(\-3 \+ 4 x\^4\)\^2 Sqrt\[2 x\^2 \+ Sqrt\[3 \+ 4 x\^4\]\]\)\} |
| nested\-nested\_11\-v | Trager | \(x\^4\-b\)/\(\(x\^4\+b\)\*Sqrt\[a\^2\*x\^4\-b\]\*Sqrt\[a\*x\^2\+Sqrt\[a\^2\*x\^4\-b\]\]\) | Problem nested\_11 \(v\) | \{"failed", "nested radicals"\} |
| nested\-nested\_11\-v | PMT | \(x\^4\-b\)/\(\(x\^4\+b\)\*Sqrt\[a\^2\*x\^4\-b\]\*Sqrt\[a\*x\^2\+Sqrt\[a\^2\*x\^4\-b\]\]\) | Problem nested\_11 \(v\) | \{"failed", "tower construction failed", \(\-b \+ x\^4\)/\(\(b \+ x\^4\) Sqrt\[\-b \+ a\^2 x\^4\] Sqrt\[a x\^2 \+ Sqrt\[\-b \+ a\^2 x\^… |
| nested\-nested\_11\-v\-inst | Trager | \(x\^4\-3\)/\(\(x\^4\+3\)\*Sqrt\[2\^2\*x\^4\-3\]\*Sqrt\[2\*x\^2\+Sqrt\[2\^2\*x\^4\-3\]\]\) | Problem nested\_11 \(v\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_11\-v\-inst | PMT | \(x\^4\-3\)/\(\(x\^4\+3\)\*Sqrt\[2\^2\*x\^4\-3\]\*Sqrt\[2\*x\^2\+Sqrt\[2\^2\*x\^4\-3\]\]\) | Problem nested\_11 \(v\); inst a=2, b=3 | \{"failed", "tower construction failed", \(\-3 \+ x\^4\)/\(\(3 \+ x\^4\) Sqrt\[\-3 \+ 4 x\^4\] Sqrt\[2 x\^2 \+ Sqrt\[\-3 \+ 4 x\^4\]\]\)… |
| nested\-nested\_12\-i | Trager | Sqrt\[a\^2\*x\^2\+b\]/\(a\*x\-Sqrt\[a\^2\*x\^2\+b\]\)\^\(1/3\) | Problem nested\_12 \(i\) | \{"failed", "nested radicals"\} |
| nested\-nested\_12\-i | PMT | Sqrt\[a\^2\*x\^2\+b\]/\(a\*x\-Sqrt\[a\^2\*x\^2\+b\]\)\^\(1/3\) | Problem nested\_12 \(i\) | \{"failed", "tower construction failed", Sqrt\[b \+ a\^2 x\^2\]/\(a x \- Sqrt\[b \+ a\^2 x\^2\]\)\^\(1/3\)\} |
| nested\-nested\_12\-i\-inst | Trager | Sqrt\[2\^2\*x\^2\+3\]/\(2\*x\-Sqrt\[2\^2\*x\^2\+3\]\)\^\(1/3\) | Problem nested\_12 \(i\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_12\-ii | Trager | \(a\*x\-Sqrt\[a\^2\*x\^2\+b\]\)\^\(1/3\)/\(x\*Sqrt\[a\^2\*x\^2\+b\]\) | Problem nested\_12 \(ii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_12\-ii | PMT | \(a\*x\-Sqrt\[a\^2\*x\^2\+b\]\)\^\(1/3\)/\(x\*Sqrt\[a\^2\*x\^2\+b\]\) | Problem nested\_12 \(ii\) | \{"failed", "tower construction failed", \(a x \- Sqrt\[b \+ a\^2 x\^2\]\)\^\(1/3\)/\(x Sqrt\[b \+ a\^2 x\^2\]\)\} |
| nested\-nested\_12\-ii\-inst | Trager | \(2\*x\-Sqrt\[2\^2\*x\^2\+3\]\)\^\(1/3\)/\(x\*Sqrt\[2\^2\*x\^2\+3\]\) | Problem nested\_12 \(ii\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_12\-iii | Trager | \(a\*x\-Sqrt\[a\^2\*x\^2\+b\]\)\^\(1/3\)/\(\(a\^2\*x\^2\+b\)\^2\*Sqrt\[a\^2\*x\^2\+b\]\) | Problem nested\_12 \(iii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_12\-iii | PMT | \(a\*x\-Sqrt\[a\^2\*x\^2\+b\]\)\^\(1/3\)/\(\(a\^2\*x\^2\+b\)\^2\*Sqrt\[a\^2\*x\^2\+b\]\) | Problem nested\_12 \(iii\) | \{"failed", "tower construction failed", \(a x \- Sqrt\[b \+ a\^2 x\^2\]\)\^\(1/3\)/\(b \+ a\^2 x\^2\)\^\(5/2\)\} |
| nested\-nested\_12\-iii\-inst | Trager | \(2\*x\-Sqrt\[2\^2\*x\^2\+3\]\)\^\(1/3\)/\(\(2\^2\*x\^2\+3\)\^2\*Sqrt\[2\^2\*x\^2\+3\]\) | Problem nested\_12 \(iii\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_12\-iv | Trager | \(a\^2\*x\^2\-b\)/\(\(a\^2\*x\^2\+b\)\*\(a\*x\-Sqrt\[a\^2\*x\^2\+b\]\)\^\(1/3\)\) | Problem nested\_12 \(iv\) | \{"failed", "nested radicals"\} |
| nested\-nested\_12\-iv | PMT | \(a\^2\*x\^2\-b\)/\(\(a\^2\*x\^2\+b\)\*\(a\*x\-Sqrt\[a\^2\*x\^2\+b\]\)\^\(1/3\)\) | Problem nested\_12 \(iv\) | \{"failed", "tower construction failed", \(\-b \+ a\^2 x\^2\)/\(\(b \+ a\^2 x\^2\) \(a x \- Sqrt\[b \+ a\^2 x\^2\]\)\^\(1/3\)\)\} |
| nested\-nested\_12\-iv\-inst | Trager | \(2\^2\*x\^2\-3\)/\(\(2\^2\*x\^2\+3\)\*\(2\*x\-Sqrt\[2\^2\*x\^2\+3\]\)\^\(1/3\)\) | Problem nested\_12 \(iv\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_13\-i | Trager | x\*Sqrt\[a\*x\^2\-b\]/Sqrt\[\(c\*x\*Sqrt\[a\*x\^2\-b\]\-a\*x\^2\+b\)\] | Problem nested\_13 \(i\) | \{"failed", "nested radicals"\} |
| nested\-nested\_13\-i | PMT | x\*Sqrt\[a\*x\^2\-b\]/Sqrt\[\(c\*x\*Sqrt\[a\*x\^2\-b\]\-a\*x\^2\+b\)\] | Problem nested\_13 \(i\) | \{"failed", "tower construction failed", \(x Sqrt\[\-b \+ a x\^2\]\)/Sqrt\[b \- a x\^2 \+ c x Sqrt\[\-b \+ a x\^2\]\]\} |
| nested\-nested\_13\-i\-inst | Trager | x\*Sqrt\[2\*x\^2\-3\]/Sqrt\[\(5\*x\*Sqrt\[2\*x\^2\-3\]\-2\*x\^2\+3\)\] | Problem nested\_13 \(i\); inst a=2, b=3, c=… | \{"failed", "nested radicals"\} |
| nested\-nested\_13\-i\-inst | PMT | x\*Sqrt\[2\*x\^2\-3\]/Sqrt\[\(5\*x\*Sqrt\[2\*x\^2\-3\]\-2\*x\^2\+3\)\] | Problem nested\_13 \(i\); inst a=2, b=3, c=… | \{"failed", "tower construction failed", \(x Sqrt\[\-3 \+ 2 x\^2\]\)/Sqrt\[3 \- 2 x\^2 \+ 5 x Sqrt\[\-3 \+ 2 x\^2\]\]\} |
| nested\-nested\_13\-ii | Trager | \(c\*x\*Sqrt\[a\*x\^2\-b\]\-a\*x\^2\+b\)\^\(3/2\)/\(x\*Sqrt\[a\*x\^2\-b\]\) | Problem nested\_13 \(ii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_13\-ii | PMT | \(c\*x\*Sqrt\[a\*x\^2\-b\]\-a\*x\^2\+b\)\^\(3/2\)/\(x\*Sqrt\[a\*x\^2\-b\]\) | Problem nested\_13 \(ii\) | \{"failed", "tower construction failed", \(b \- a x\^2 \+ c x Sqrt\[\-b \+ a x\^2\]\)\^\(3/2\)/\(x Sqrt\[\-b \+ a x\^2\]\)\} |
| nested\-nested\_13\-ii\-inst | Trager | \(5\*x\*Sqrt\[2\*x\^2\-3\]\-2\*x\^2\+3\)\^\(3/2\)/\(x\*Sqrt\[2\*x\^2\-3\]\) | Problem nested\_13 \(ii\); inst a=2, b=3, c… | \{"failed", "nested radicals"\} |
| nested\-nested\_13\-ii\-inst | PMT | \(5\*x\*Sqrt\[2\*x\^2\-3\]\-2\*x\^2\+3\)\^\(3/2\)/\(x\*Sqrt\[2\*x\^2\-3\]\) | Problem nested\_13 \(ii\); inst a=2, b=3, c… | \{"failed", "tower construction failed", \(3 \- 2 x\^2 \+ 5 x Sqrt\[\-3 \+ 2 x\^2\]\)\^\(3/2\)/\(x Sqrt\[\-3 \+ 2 x\^2\]\)\} |
| nested\-nested\_13\-iii | Trager | Sqrt\[\(c\*x\*Sqrt\[a\*x\^2\-b\]\-a\*x\^2\+b\)\]/\(c\*x\*Sqrt\[a\*x\^2\-b\]\-a\*x\^2\) | Problem nested\_13 \(iii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_13\-iii | PMT | Sqrt\[\(c\*x\*Sqrt\[a\*x\^2\-b\]\-a\*x\^2\+b\)\]/\(c\*x\*Sqrt\[a\*x\^2\-b\]\-a\*x\^2\) | Problem nested\_13 \(iii\) | \{"failed", "tower construction failed", Sqrt\[b \- a x\^2 \+ c x Sqrt\[\-b \+ a x\^2\]\]/\(\-a x\^2 \+ c x Sqrt\[\-b \+ a x\^2\]\)… |
| nested\-nested\_13\-iii\-inst | Trager | Sqrt\[\(5\*x\*Sqrt\[2\*x\^2\-3\]\-2\*x\^2\+3\)\]/\(5\*x\*Sqrt\[2\*x\^2\-3\]\-2\*x\^2\) | Problem nested\_13 \(iii\); inst a=2, b=3, … | \{"failed", "nested radicals"\} |
| nested\-nested\_13\-iii\-inst | PMT | Sqrt\[\(5\*x\*Sqrt\[2\*x\^2\-3\]\-2\*x\^2\+3\)\]/\(5\*x\*Sqrt\[2\*x\^2\-3\]\-2\*x\^2\) | Problem nested\_13 \(iii\); inst a=2, b=3, … | \{"failed", "tower construction failed", Sqrt\[3 \- 2 x\^2 \+ 5 x Sqrt\[\-3 \+ 2 x\^2\]\]/\(\-2 x\^2 \+ 5 x Sqrt\[\-3 \+ 2 x\^2\]\)… |
| nested\-nested\_13\-iv | Trager | \(c\*x\*Sqrt\[a\*x\^2\-b\]\+b\)/\(c\*x\*Sqrt\[a\*x\^2\-b\]\-a\*x\^2\+b\)\^\(5/2\) | Problem nested\_13 \(iv\) | \{"failed", "nested radicals"\} |
| nested\-nested\_13\-iv | PMT | \(c\*x\*Sqrt\[a\*x\^2\-b\]\+b\)/\(c\*x\*Sqrt\[a\*x\^2\-b\]\-a\*x\^2\+b\)\^\(5/2\) | Problem nested\_13 \(iv\) | \{"failed", "tower construction failed", \(b \+ c x Sqrt\[\-b \+ a x\^2\]\)/\(b \- a x\^2 \+ c x Sqrt\[\-b \+ a x\^2\]\)\^\(5/2\)\} |
| nested\-nested\_13\-iv\-inst | Trager | \(5\*x\*Sqrt\[2\*x\^2\-3\]\+3\)/\(5\*x\*Sqrt\[2\*x\^2\-3\]\-2\*x\^2\+3\)\^\(5/2\) | Problem nested\_13 \(iv\); inst a=2, b=3, c… | \{"failed", "nested radicals"\} |
| nested\-nested\_13\-iv\-inst | PMT | \(5\*x\*Sqrt\[2\*x\^2\-3\]\+3\)/\(5\*x\*Sqrt\[2\*x\^2\-3\]\-2\*x\^2\+3\)\^\(5/2\) | Problem nested\_13 \(iv\); inst a=2, b=3, c… | \{"failed", "tower construction failed", \(3 \+ 5 x Sqrt\[\-3 \+ 2 x\^2\]\)/\(3 \- 2 x\^2 \+ 5 x Sqrt\[\-3 \+ 2 x\^2\]\)\^\(5/2\)\} |
| nested\-nested\_13\-v | Trager | \(c\*x\*Sqrt\[a\*x\^2\-b\]\-a\*x\^2\+b\)\^\(3/2\) | Problem nested\_13 \(v\) | \{"failed", "nested radicals"\} |
| nested\-nested\_13\-v | PMT | \(c\*x\*Sqrt\[a\*x\^2\-b\]\-a\*x\^2\+b\)\^\(3/2\) | Problem nested\_13 \(v\) | \{"failed", "tower construction failed", \(b \- a x\^2 \+ c x Sqrt\[\-b \+ a x\^2\]\)\^\(3/2\)\} |
| nested\-nested\_13\-v\-inst | Trager | \(5\*x\*Sqrt\[2\*x\^2\-3\]\-2\*x\^2\+3\)\^\(3/2\) | Problem nested\_13 \(v\); inst a=2, b=3, c=… | \{"failed", "nested radicals"\} |
| nested\-nested\_13\-v\-inst | PMT | \(5\*x\*Sqrt\[2\*x\^2\-3\]\-2\*x\^2\+3\)\^\(3/2\) | Problem nested\_13 \(v\); inst a=2, b=3, c=… | \{"failed", "tower construction failed", \(3 \- 2 x\^2 \+ 5 x Sqrt\[\-3 \+ 2 x\^2\]\)\^\(3/2\)\} |
| nested\-nested\_14\-i | Trager | Sqrt\[Sqrt\[a\*x\^2\+b\^2\]\+b\] | Problem nested\_14 \(i\) | \{"failed", "nested radicals"\} |
| nested\-nested\_14\-i | PMT | Sqrt\[Sqrt\[a\*x\^2\+b\^2\]\+b\] | Problem nested\_14 \(i\) | \{"failed", "tower construction failed", Sqrt\[b \+ Sqrt\[b\^2 \+ a x\^2\]\]\} |
| nested\-nested\_14\-i\-inst | Trager | Sqrt\[Sqrt\[2\*x\^2\+3\^2\]\+3\] | Problem nested\_14 \(i\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_14\-ii | Trager | Sqrt\[Sqrt\[a\*x\^2\+b\^2\]\+b\]/Sqrt\[a\*x\^2\+b\^2\] | Problem nested\_14 \(ii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_14\-ii | PMT | Sqrt\[Sqrt\[a\*x\^2\+b\^2\]\+b\]/Sqrt\[a\*x\^2\+b\^2\] | Problem nested\_14 \(ii\) | \{"failed", "tower construction failed", Sqrt\[b \+ Sqrt\[b\^2 \+ a x\^2\]\]/Sqrt\[b\^2 \+ a x\^2\]\} |
| nested\-nested\_14\-ii\-inst | Trager | Sqrt\[Sqrt\[2\*x\^2\+3\^2\]\+3\]/Sqrt\[2\*x\^2\+3\^2\] | Problem nested\_14 \(ii\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_14\-iii | Trager | Sqrt\[Sqrt\[a\*x\^2\+b\^2\]\+b\]/\(x\^2\*Sqrt\[a\*x\^2\+b\^2\]\) | Problem nested\_14 \(iii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_14\-iii | PMT | Sqrt\[Sqrt\[a\*x\^2\+b\^2\]\+b\]/\(x\^2\*Sqrt\[a\*x\^2\+b\^2\]\) | Problem nested\_14 \(iii\) | \{"failed", "tower construction failed", Sqrt\[b \+ Sqrt\[b\^2 \+ a x\^2\]\]/\(x\^2 Sqrt\[b\^2 \+ a x\^2\]\)\} |
| nested\-nested\_14\-iii\-inst | Trager | Sqrt\[Sqrt\[2\*x\^2\+3\^2\]\+3\]/\(x\^2\*Sqrt\[2\*x\^2\+3\^2\]\) | Problem nested\_14 \(iii\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_14\-iii\-inst | PMT | Sqrt\[Sqrt\[2\*x\^2\+3\^2\]\+3\]/\(x\^2\*Sqrt\[2\*x\^2\+3\^2\]\) | Problem nested\_14 \(iii\); inst a=2, b=3 | \{"failed", "no solution within bounds", \{1\}\} |
| nested\-nested\_14\-iv | Trager | Sqrt\[Sqrt\[a\*x\^2\+b\^2\]\+b\]/\(\(a\*x\^2\-b\)\*Sqrt\[a\*x\^2\+b\^2\]\) | Problem nested\_14 \(iv\) | \{"failed", "nested radicals"\} |
| nested\-nested\_14\-iv | PMT | Sqrt\[Sqrt\[a\*x\^2\+b\^2\]\+b\]/\(\(a\*x\^2\-b\)\*Sqrt\[a\*x\^2\+b\^2\]\) | Problem nested\_14 \(iv\) | \{"failed", "tower construction failed", Sqrt\[b \+ Sqrt\[b\^2 \+ a x\^2\]\]/\(\(\-b \+ a x\^2\) Sqrt\[b\^2 \+ a x\^2\]\)\} |
| nested\-nested\_14\-iv\-inst | Trager | Sqrt\[Sqrt\[2\*x\^2\+3\^2\]\+3\]/\(\(2\*x\^2\-3\)\*Sqrt\[2\*x\^2\+3\^2\]\) | Problem nested\_14 \(iv\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_14\-v | Trager | \(a\*x\^2\+b\)\^2\*Sqrt\[a\*x\^2\+b\^2\]/\(\(a\*x\^2\-b\)\^2\*Sqrt\[Sqrt\[a\*x\^2\+b\^2\]\+b\]\) | Problem nested\_14 \(v\) | \{"failed", "nested radicals"\} |
| nested\-nested\_14\-v | PMT | \(a\*x\^2\+b\)\^2\*Sqrt\[a\*x\^2\+b\^2\]/\(\(a\*x\^2\-b\)\^2\*Sqrt\[Sqrt\[a\*x\^2\+b\^2\]\+b\]\) | Problem nested\_14 \(v\) | \{"failed", "tower construction failed", \(\(b \+ a x\^2\)\^2 Sqrt\[b\^2 \+ a x\^2\]\)/\(\(\-b \+ a x\^2\)\^2 Sqrt\[b \+ Sqrt\[b\^2 \+ … |
| nested\-nested\_14\-v\-inst | Trager | \(2\*x\^2\+3\)\^2\*Sqrt\[2\*x\^2\+3\^2\]/\(\(2\*x\^2\-3\)\^2\*Sqrt\[Sqrt\[2\*x\^2\+3\^2\]\+3\]\) | Problem nested\_14 \(v\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_15\-i | Trager | Sqrt\[Sqrt\[a\^2\*x\^2\+b\]\+a\*x\] | Problem nested\_15 \(i\) | \{"failed", "nested radicals"\} |
| nested\-nested\_15\-i | PMT | Sqrt\[Sqrt\[a\^2\*x\^2\+b\]\+a\*x\] | Problem nested\_15 \(i\) | \{"failed", "tower construction failed", Sqrt\[a x \+ Sqrt\[b \+ a\^2 x\^2\]\]\} |
| nested\-nested\_15\-i\-inst | Trager | Sqrt\[Sqrt\[2\^2\*x\^2\+3\]\+2\*x\] | Problem nested\_15 \(i\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_15\-ii | Trager | x\^2\*Sqrt\[Sqrt\[a\^2\*x\^2\+b\]\+a\*x\] | Problem nested\_15 \(ii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_15\-ii | PMT | x\^2\*Sqrt\[Sqrt\[a\^2\*x\^2\+b\]\+a\*x\] | Problem nested\_15 \(ii\) | \{"failed", "tower construction failed", x\^2 Sqrt\[a x \+ Sqrt\[b \+ a\^2 x\^2\]\]\} |
| nested\-nested\_15\-ii\-inst | Trager | x\^2\*Sqrt\[Sqrt\[2\^2\*x\^2\+3\]\+2\*x\] | Problem nested\_15 \(ii\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_15\-iii | Trager | 1/\(x\^2\*Sqrt\[Sqrt\[a\^2\*x\^2\+b\]\+a\*x\]\) | Problem nested\_15 \(iii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_15\-iii | PMT | 1/\(x\^2\*Sqrt\[Sqrt\[a\^2\*x\^2\+b\]\+a\*x\]\) | Problem nested\_15 \(iii\) | \{"failed", "tower construction failed", 1/\(x\^2 Sqrt\[a x \+ Sqrt\[b \+ a\^2 x\^2\]\]\)\} |
| nested\-nested\_15\-iii\-inst | Trager | 1/\(x\^2\*Sqrt\[Sqrt\[2\^2\*x\^2\+3\]\+2\*x\]\) | Problem nested\_15 \(iii\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_15\-iv | Trager | \(a\^2\*x\^2\+b\)\^\(3/2\)/\(x\*\(Sqrt\[a\^2\*x\^2\+b\]\+a\*x\)\^\(3/2\)\) | Problem nested\_15 \(iv\) | \{"failed", "nested radicals"\} |
| nested\-nested\_15\-iv | PMT | \(a\^2\*x\^2\+b\)\^\(3/2\)/\(x\*\(Sqrt\[a\^2\*x\^2\+b\]\+a\*x\)\^\(3/2\)\) | Problem nested\_15 \(iv\) | \{"failed", "tower construction failed", \(b \+ a\^2 x\^2\)\^\(3/2\)/\(x \(a x \+ Sqrt\[b \+ a\^2 x\^2\]\)\^\(3/2\)\)\} |
| nested\-nested\_15\-iv\-inst | Trager | \(2\^2\*x\^2\+3\)\^\(3/2\)/\(x\*\(Sqrt\[2\^2\*x\^2\+3\]\+2\*x\)\^\(3/2\)\) | Problem nested\_15 \(iv\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_15\-v | Trager | x\^2/\(\(a\^2\*x\^2\+b\)\^\(3/2\)\*\(Sqrt\[a\^2\*x\^2\+b\]\+a\*x\)\^\(3/2\)\) | Problem nested\_15 \(v\) | \{"failed", "nested radicals"\} |
| nested\-nested\_15\-v | PMT | x\^2/\(\(a\^2\*x\^2\+b\)\^\(3/2\)\*\(Sqrt\[a\^2\*x\^2\+b\]\+a\*x\)\^\(3/2\)\) | Problem nested\_15 \(v\) | \{"failed", "tower construction failed", x\^2/\(\(b \+ a\^2 x\^2\)\^\(3/2\) \(a x \+ Sqrt\[b \+ a\^2 x\^2\]\)\^\(3/2\)\)\} |
| nested\-nested\_15\-v\-inst | Trager | x\^2/\(\(2\^2\*x\^2\+3\)\^\(3/2\)\*\(Sqrt\[2\^2\*x\^2\+3\]\+2\*x\)\^\(3/2\)\) | Problem nested\_15 \(v\); inst a=2, b=3 | \{"failed", "nested radicals"\} |
| nested\-nested\_16\-i | Trager | Sqrt\[Sqrt\[Sqrt\[a\^2\*x\^2\+b\]\+a\*x\]\+c\] | Problem nested\_16 \(i\) | \{"failed", "nested radicals"\} |
| nested\-nested\_16\-i | PMT | Sqrt\[Sqrt\[Sqrt\[a\^2\*x\^2\+b\]\+a\*x\]\+c\] | Problem nested\_16 \(i\) | \{"failed", "tower construction failed", Sqrt\[c \+ Sqrt\[a x \+ Sqrt\[b \+ a\^2 x\^2\]\]\]\} |
| nested\-nested\_16\-i\-inst | Trager | Sqrt\[Sqrt\[Sqrt\[2\^2\*x\^2\+3\]\+2\*x\]\+5\] | Problem nested\_16 \(i\); inst a=2, b=3, c=… | \{"failed", "nested radicals"\} |
| nested\-nested\_16\-ii | Trager | Sqrt\[Sqrt\[Sqrt\[a\^2\*x\^2\+b\]\+a\*x\]\+c\]/\(Sqrt\[a\^2\*x\^2\+b\]\*Sqrt\[Sqrt\[a\^2\*x\^2\+b\]\+a\*x\]\) | Problem nested\_16 \(ii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_16\-ii | PMT | Sqrt\[Sqrt\[Sqrt\[a\^2\*x\^2\+b\]\+a\*x\]\+c\]/\(Sqrt\[a\^2\*x\^2\+b\]\*Sqrt\[Sqrt\[a\^2\*x\^2\+b\]\+a\*x\]\) | Problem nested\_16 \(ii\) | \{"failed", "tower construction failed", Sqrt\[c \+ Sqrt\[a x \+ Sqrt\[b \+ a\^2 x\^2\]\]\]/\(Sqrt\[b \+ a\^2 x\^2\] Sqrt\[a x \+ … |
| nested\-nested\_16\-ii\-inst | Trager | Sqrt\[Sqrt\[Sqrt\[2\^2\*x\^2\+3\]\+2\*x\]\+5\]/\(Sqrt\[2\^2\*x\^2\+3\]\*Sqrt\[Sqrt\[2\^2\*x\^2\+3\]\+2\*x\]\) | Problem nested\_16 \(ii\); inst a=2, b=3, c… | \{"failed", "nested radicals"\} |
| nested\-nested\_16\-iii | Trager | Sqrt\[Sqrt\[Sqrt\[a\^2\*x\^2\+b\]\+a\*x\]\+c\]/Sqrt\[Sqrt\[a\^2\*x\^2\+b\]\+a\*x\] | Problem nested\_16 \(iii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_16\-iii | PMT | Sqrt\[Sqrt\[Sqrt\[a\^2\*x\^2\+b\]\+a\*x\]\+c\]/Sqrt\[Sqrt\[a\^2\*x\^2\+b\]\+a\*x\] | Problem nested\_16 \(iii\) | \{"failed", "tower construction failed", Sqrt\[c \+ Sqrt\[a x \+ Sqrt\[b \+ a\^2 x\^2\]\]\]/Sqrt\[a x \+ Sqrt\[b \+ a\^2 x\^2\]\]\} |
| nested\-nested\_16\-iii\-inst | Trager | Sqrt\[Sqrt\[Sqrt\[2\^2\*x\^2\+3\]\+2\*x\]\+5\]/Sqrt\[Sqrt\[2\^2\*x\^2\+3\]\+2\*x\] | Problem nested\_16 \(iii\); inst a=2, b=3, … | \{"failed", "nested radicals"\} |
| nested\-nested\_17\-i | Trager | \(\(x\^2\+1\)\^\(3/2\)\+x\^3\)\^\(1/3\)/\(x\^2\*\(x\^2\+1\)\) | Problem nested\_17 \(i\) | \{"failed", "nested radicals"\} |
| nested\-nested\_17\-ii | Trager | \(a\*\(c\*x\^2\+d\)\^\(3/2\)\+b\*x\^3\)\^\(1/3\)/\(x\^2\*\(c\*x\^2\+d\)\) | Problem nested\_17 \(ii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_17\-ii | PMT | \(a\*\(c\*x\^2\+d\)\^\(3/2\)\+b\*x\^3\)\^\(1/3\)/\(x\^2\*\(c\*x\^2\+d\)\) | Problem nested\_17 \(ii\) | \{"failed", "tower construction failed", \(b x\^3 \+ a \(d \+ c x\^2\)\^\(3/2\)\)\^\(1/3\)/\(x\^2 \(d \+ c x\^2\)\)\} |
| nested\-nested\_17\-ii\-inst | Trager | \(2\*\(5\*x\^2\+7\)\^\(3/2\)\+3\*x\^3\)\^\(1/3\)/\(x\^2\*\(5\*x\^2\+7\)\) | Problem nested\_17 \(ii\); inst a=2, b=3, c… | \{"failed", "nested radicals"\} |
| nested\-nested\_17\-ii\-inst | PMT | \(2\*\(5\*x\^2\+7\)\^\(3/2\)\+3\*x\^3\)\^\(1/3\)/\(x\^2\*\(5\*x\^2\+7\)\) | Problem nested\_17 \(ii\); inst a=2, b=3, c… | \{"failed", "tower construction failed", \(3 x\^3 \+ 2 \(7 \+ 5 x\^2\)\^\(3/2\)\)\^\(1/3\)/\(x\^2 \(7 \+ 5 x\^2\)\)\} |
| nested\-nested\_17\-iii | Trager | \(2\*x\^2\+1\)\*\(x\^3\*\(x\^2\+1\)\^\(3/2\)\+1\)\^\(1/3\)/\(x\*\(x\^2\+1\)\) | Problem nested\_17 \(iii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_17\-iv | Trager | \(2\*c\*x\^2\+d\)\*\(a\*x\^3\*\(c\*x\^2\+d\)\^\(3/2\)\+b\)\^\(1/3\)/\(x\*\(c\*x\^2\+d\)\) | Problem nested\_17 \(iv\) | \{"failed", "nested radicals"\} |
| nested\-nested\_17\-iv | PMT | \(2\*c\*x\^2\+d\)\*\(a\*x\^3\*\(c\*x\^2\+d\)\^\(3/2\)\+b\)\^\(1/3\)/\(x\*\(c\*x\^2\+d\)\) | Problem nested\_17 \(iv\) | \{"failed", "tower construction failed", \(\(d \+ 2 c x\^2\) \(b \+ a x\^3 \(d \+ c x\^2\)\^\(3/2\)\)\^\(1/3\)\)/\(x \(d \+ c x\^2\)\)\} |
| nested\-nested\_17\-iv\-inst | Trager | \(2\*5\*x\^2\+7\)\*\(2\*x\^3\*\(5\*x\^2\+7\)\^\(3/2\)\+3\)\^\(1/3\)/\(x\*\(5\*x\^2\+7\)\) | Problem nested\_17 \(iv\); inst a=2, b=3, c… | \{"failed", "nested radicals"\} |
| nested\-nested\_17\-iv\-inst | PMT | \(2\*5\*x\^2\+7\)\*\(2\*x\^3\*\(5\*x\^2\+7\)\^\(3/2\)\+3\)\^\(1/3\)/\(x\*\(5\*x\^2\+7\)\) | Problem nested\_17 \(iv\); inst a=2, b=3, c… | \{"failed", "tower construction failed", \(\(7 \+ 10 x\^2\) \(3 \+ 2 x\^3 \(7 \+ 5 x\^2\)\^\(3/2\)\)\^\(1/3\)\)/\(x \(7 \+ 5 x\^2\)\)\} |
| nested\-nested\_17\-v | Trager | x\^2\*\(5\*x\^3\+2\)\*Sqrt\[x\^3\+1\]/\(x\^3\*\(x\^3\+1\)\^\(3/2\)\+1\)\^\(1/3\) | Problem nested\_17 \(v\) | \{"failed", "nested radicals"\} |
| nested\-nested\_17\-v | PMT | x\^2\*\(5\*x\^3\+2\)\*Sqrt\[x\^3\+1\]/\(x\^3\*\(x\^3\+1\)\^\(3/2\)\+1\)\^\(1/3\) | Problem nested\_17 \(v\) | \{"failed", "tower construction failed", \(x\^2 Sqrt\[1 \+ x\^3\] \(2 \+ 5 x\^3\)\)/\(1 \+ x\^3 \(1 \+ x\^3\)\^\(3/2\)\)\^\(1/3\)\} |
| nested\-nested\_17\-vi | Trager | x\^2\*\(5\*c\*x\^3\+2\*d\)\*Sqrt\[c\*x\^3\+d\]/\(a\*x\^3\*\(c\*x\^3\+d\)\^\(3/2\)\+b\)\^\(1/3\) | Problem nested\_17 \(vi\) | \{"failed", "nested radicals"\} |
| nested\-nested\_17\-vi | PMT | x\^2\*\(5\*c\*x\^3\+2\*d\)\*Sqrt\[c\*x\^3\+d\]/\(a\*x\^3\*\(c\*x\^3\+d\)\^\(3/2\)\+b\)\^\(1/3\) | Problem nested\_17 \(vi\) | \{"failed", "tower construction failed", \(x\^2 Sqrt\[d \+ c x\^3\] \(2 d \+ 5 c x\^3\)\)/\(b \+ a x\^3 \(d \+ c x\^3\)\^\(3/2\)\)\^\(1… |
| nested\-nested\_17\-vi\-inst | Trager | x\^2\*\(5\*5\*x\^3\+2\*7\)\*Sqrt\[5\*x\^3\+7\]/\(2\*x\^3\*\(5\*x\^3\+7\)\^\(3/2\)\+3\)\^\(1/3\) | Problem nested\_17 \(vi\); inst a=2, b=3, c… | \{"failed", "nested radicals"\} |
| nested\-nested\_17\-vi\-inst | PMT | x\^2\*\(5\*5\*x\^3\+2\*7\)\*Sqrt\[5\*x\^3\+7\]/\(2\*x\^3\*\(5\*x\^3\+7\)\^\(3/2\)\+3\)\^\(1/3\) | Problem nested\_17 \(vi\); inst a=2, b=3, c… | \{"failed", "tower construction failed", \(x\^2 Sqrt\[7 \+ 5 x\^3\] \(14 \+ 25 x\^3\)\)/\(3 \+ 2 x\^3 \(7 \+ 5 x\^3\)\^\(3/2\)\)\^\(1/3… |
| nested\-nested\_19\-i | Trager | 1/Sqrt\[x\^4\+x\^2\*Sqrt\[x\^4\+4\]\] | Problem nested\_19 \(i\) | \{"failed", "nested radicals"\} |
| nested\-nested\_19\-i | PMT | 1/Sqrt\[x\^4\+x\^2\*Sqrt\[x\^4\+4\]\] | Problem nested\_19 \(i\) | \{"failed", "tower construction failed", 1/Sqrt\[x\^4 \+ x\^2 Sqrt\[4 \+ x\^4\]\]\} |
| nested\-nested\_19\-ii | Trager | Sqrt\[1/x\^4\-Sqrt\[4\*x\^4\+1\]/x\^4\] | Problem nested\_19 \(ii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_19\-ii | PMT | Sqrt\[1/x\^4\-Sqrt\[4\*x\^4\+1\]/x\^4\] | Problem nested\_19 \(ii\) | \{"failed", "tower construction failed", Sqrt\[1/x\^4 \- Sqrt\[1 \+ 4 x\^4\]/x\^4\]\} |
| nested\-nested\_19\-iii | Trager | \(Sqrt\[2\]\*x\^3\-Sqrt\[2\*x\^6\-1\]\)\^\(1/3\) | Problem nested\_19 \(iii\) | \{"failed", "algebraic constant in the integrand \(coefficients must be rational\)"\} |
| nested\-nested\_19\-iii | PMT | \(Sqrt\[2\]\*x\^3\-Sqrt\[2\*x\^6\-1\]\)\^\(1/3\) | Problem nested\_19 \(iii\) | \{"failed", "tower construction failed", \(Sqrt\[2\] x\^3 \- Sqrt\[\-1 \+ 2 x\^6\]\)\^\(1/3\)\} |
| nested\-nested\_19\-iv | Trager | \(Sqrt\[2\]/x\^3\-Sqrt\[2\-x\^6\]/x\^3\)\^\(1/3\) | Problem nested\_19 \(iv\) | \{"failed", "algebraic constant in the integrand \(coefficients must be rational\)"\} |
| nested\-nested\_19\-iv | PMT | \(Sqrt\[2\]/x\^3\-Sqrt\[2\-x\^6\]/x\^3\)\^\(1/3\) | Problem nested\_19 \(iv\) | \{"failed", "tower construction failed", \(Sqrt\[2\]/x\^3 \- Sqrt\[2 \- x\^6\]/x\^3\)\^\(1/3\)\} |
| nested\-nested\_19\-v | Trager | \(x\^3\+Sqrt\[x\^6\-1\]\)\^\(1/3\) | Problem nested\_19 \(v\) | \{"failed", "nested radicals"\} |
| nested\-nested\_19\-v | PMT | \(x\^3\+Sqrt\[x\^6\-1\]\)\^\(1/3\) | Problem nested\_19 \(v\) | \{"failed", "tower construction failed", \(x\^3 \+ Sqrt\[\-1 \+ x\^6\]\)\^\(1/3\)\} |
| nested\-nested\_19\-vi | Trager | \(1/x\^3\+Sqrt\[1/x\^6\-4\]\)\^\(1/3\) | Problem nested\_19 \(vi\) | \{"failed", "nested radicals"\} |
| nested\-nested\_19\-vi | PMT | \(1/x\^3\+Sqrt\[1/x\^6\-4\]\)\^\(1/3\) | Problem nested\_19 \(vi\) | \{"failed", "tower construction failed", \(Sqrt\[\-4 \+ 1/x\^6\] \+ 1/x\^3\)\^\(1/3\)\} |
| nested\-nested\_19\-vii | Trager | \(1\+x\^3\+Sqrt\[x\^6\+2\*x\^3\]\)\^\(1/3\) | Problem nested\_19 \(vii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_19\-vii | PMT | \(1\+x\^3\+Sqrt\[x\^6\+2\*x\^3\]\)\^\(1/3\) | Problem nested\_19 \(vii\) | \{"failed", "tower construction failed", \(1 \+ x\^3 \+ Sqrt\[2 x\^3 \+ x\^6\]\)\^\(1/3\)\} |
| nested\-nested\_19\-viii | Trager | \(1\+1/x\^3\+Sqrt\[2\*x\^3\+1\]/x\^3\)\^\(1/3\) | Problem nested\_19 \(viii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_19\-viii | PMT | \(1\+1/x\^3\+Sqrt\[2\*x\^3\+1\]/x\^3\)\^\(1/3\) | Problem nested\_19 \(viii\) | \{"failed", "tower construction failed", \(1 \+ 1/x\^3 \+ Sqrt\[1 \+ 2 x\^3\]/x\^3\)\^\(1/3\)\} |
| nested\-nested\_19\-ix | Trager | \(x\^3\-Sqrt\[x\^6\+4\]\)\^\(1/3\) | Problem nested\_19 \(ix\) | \{"failed", "nested radicals"\} |
| nested\-nested\_19\-ix | PMT | \(x\^3\-Sqrt\[x\^6\+4\]\)\^\(1/3\) | Problem nested\_19 \(ix\) | \{"failed", "tower construction failed", \(x\^3 \- Sqrt\[4 \+ x\^6\]\)\^\(1/3\)\} |
| nested\-nested\_19\-x | Trager | \(1/x\^3\-Sqrt\[4\*x\^6\+1\]/x\^3\)\^\(1/3\) | Problem nested\_19 \(x\) | \{"failed", "nested radicals"\} |
| nested\-nested\_19\-x | PMT | \(1/x\^3\-Sqrt\[4\*x\^6\+1\]/x\^3\)\^\(1/3\) | Problem nested\_19 \(x\) | \{"failed", "tower construction failed", \(1/x\^3 \- Sqrt\[1 \+ 4 x\^6\]/x\^3\)\^\(1/3\)\} |
| nested\-nested\_19\-xi | Trager | \(2\+1/x\^2\-Sqrt\[8\*x\^2\+1\]/x\^2\)\^\(1/4\) | Problem nested\_19 \(xi\) | \{"failed", "nested radicals"\} |
| nested\-nested\_19\-xii | Trager | \(1\-Sqrt\[x\^8\+8\*x\^4\]/x\^4\)\^\(1/4\) | Problem nested\_19 \(xii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_19\-xii | PMT | \(1\-Sqrt\[x\^8\+8\*x\^4\]/x\^4\)\^\(1/4\) | Problem nested\_19 \(xii\) | \{"failed", "tower construction failed", \(1 \- Sqrt\[8 x\^4 \+ x\^8\]/x\^4\)\^\(1/4\)\} |
| nested\-nested\_19\-xiii | Trager | \(1\+Sqrt\[x\^8\+1\]/x\^4\)\^\(1/4\) | Problem nested\_19 \(xiii\) | \{"failed", "nested radicals"\} |
| nested\-nested\_19\-xiii | PMT | \(1\+Sqrt\[x\^8\+1\]/x\^4\)\^\(1/4\) | Problem nested\_19 \(xiii\) | \{"failed", "tower construction failed", \(1 \+ Sqrt\[1 \+ x\^8\]/x\^4\)\^\(1/4\)\} |
| multrad\-intro\-1 | Trager | \(x\^4\-12\)/\(\(x\^4\-8\*x\^2\+12\)\*Sqrt\[x\^2\-4\]\*Sqrt\[x\^2\-3\]\) | Section intro example \(solved, \-ArcTanh\) | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_0\-i | Trager | \(x\^4\-a\*b\)/\(\(d\*x\^4\-\(a\*d\+b\*d\-c\)\*x\^2\+a\*b\*d\)\*Sqrt\[x\^2\-a\]\*Sqrt\[x\^2\-b\]\) | Problem multiple\_radicals\_0 \(i\) | \{"failed", "radicand is not a rational function over Q"\} |
| multrad\-multiple\_radicals\_0\-i\-inst | Trager | \(x\^4\-2\*3\)/\(\(7\*x\^4\-\(2\*7\+3\*7\-5\)\*x\^2\+2\*3\*7\)\*Sqrt\[x\^2\-2\]\*Sqrt\[x\^2\-3\]\) | Problem multiple\_radicals\_0 \(i\); inst a=… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_0\-ii | Trager | \(x\^4\-a\*b\)/\(\(x\^4\+a\*b\)\*Sqrt\[x\^2\-a\]\*Sqrt\[x\^2\-b\]\) | Problem multiple\_radicals\_0 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| multrad\-multiple\_radicals\_0\-ii\-inst | Trager | \(x\^4\-2\*3\)/\(\(x\^4\+2\*3\)\*Sqrt\[x\^2\-2\]\*Sqrt\[x\^2\-3\]\) | Problem multiple\_radicals\_0 \(ii\); inst a… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_0\-iii | Trager | \(a\*x\^4\+b\)\^\(1/4\)/\(b\^2\-a\^2\*x\^8\)\^\(1/4\) | Problem multiple\_radicals\_0 \(iii\) | \{"failed", "radicand is not a rational function over Q"\} |
| multrad\-multiple\_radicals\_0\-iii | PMT | \(a\*x\^4\+b\)\^\(1/4\)/\(b\^2\-a\^2\*x\^8\)\^\(1/4\) | Problem multiple\_radicals\_0 \(iii\) | \{"failed", "tower construction failed", \(b \+ a x\^4\)\^\(1/4\)/\(b\^2 \- a\^2 x\^8\)\^\(1/4\)\} |
| multrad\-multiple\_radicals\_0\-iii\-inst | Trager | \(2\*x\^4\+3\)\^\(1/4\)/\(3\^2\-2\^2\*x\^8\)\^\(1/4\) | Problem multiple\_radicals\_0 \(iii\); inst … | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_0\-iii\-inst | PMT | \(2\*x\^4\+3\)\^\(1/4\)/\(3\^2\-2\^2\*x\^8\)\^\(1/4\) | Problem multiple\_radicals\_0 \(iii\); inst … | \{"failed", "tower construction failed", \(3 \+ 2 x\^4\)\^\(1/4\)/\(9 \- 4 x\^8\)\^\(1/4\)\} |
| multrad\-multiple\_radicals\_0\-iv | Trager | \(x\^4\-a\*b\)\*Sqrt\[x\^2\-a\]\*Sqrt\[x\^2\-b\]/\(x\^2\*\(d\*x\^4\+\(c\-a\*d\-b\*d\)\*x\^2\+a\*b\*d\)\) | Problem multiple\_radicals\_0 \(iv\) | \{"failed", "radicand is not a rational function over Q"\} |
| multrad\-multiple\_radicals\_0\-iv\-inst | Trager | \(x\^4\-2\*3\)\*Sqrt\[x\^2\-2\]\*Sqrt\[x\^2\-3\]/\(x\^2\*\(7\*x\^4\+\(5\-2\*7\-3\*7\)\*x\^2\+2\*3\*7\)\) | Problem multiple\_radicals\_0 \(iv\); inst a… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_0\-v | Trager | x\^2\*\(x\^4\-2\*b\*x\^2\+a\*b\)\*Sqrt\[x\^2\-a\]/\(\(d\*x\^8\-2\*a\*d\*x\^6\+\(a\^2\*d\+c\)\*x\^4\-2\*b\*c\*x\^2\+b\^2\*c\)\*Sqrt\[x\^… | Problem multiple\_radicals\_0 \(v\) | \{"failed", "radicand is not a rational function over Q"\} |
| multrad\-multiple\_radicals\_0\-v\-inst | Trager | x\^2\*\(x\^4\-2\*3\*x\^2\+2\*3\)\*Sqrt\[x\^2\-2\]/\(\(7\*x\^8\-2\*2\*7\*x\^6\+\(2\^2\*7\+5\)\*x\^4\-2\*3\*5\*x\^2\+3\^2\*5\)\*Sqrt\[x\^… | Problem multiple\_radicals\_0 \(v\); inst a=… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_0\-v\-inst | PMT | x\^2\*\(x\^4\-2\*3\*x\^2\+2\*3\)\*Sqrt\[x\^2\-2\]/\(\(7\*x\^8\-2\*2\*7\*x\^6\+\(2\^2\*7\+5\)\*x\^4\-2\*3\*5\*x\^2\+3\^2\*5\)\*Sqrt\[x\^… | Problem multiple\_radicals\_0 \(v\); inst a=… | \{"failed", "residue class lost: no places over the prime for a root of the residue polynomial", 45927 \+ 40824 … |
| multrad\-multiple\_radicals\_0\-vi | Trager | \(x\^4\+2\*a\*x\^2\+a\^2\)\*\(x\^6\+3\*a\*x\^4\-a\*\(2\*a\+3\*b\)\*x\^2\+a\^2\*b\)/\(\(c\*x\^12\-2\*c\*x\^10\+\(a\^2\*c\+4\*a\*b\*c\+b\^2… | Problem multiple\_radicals\_0 \(vi\) | \{"failed", "radicand is not a rational function over Q"\} |
| multrad\-multiple\_radicals\_0\-vi\-inst | Trager | \(x\^4\+2\*2\*x\^2\+2\^2\)\*\(x\^6\+3\*2\*x\^4\-2\*\(2\*2\+3\*3\)\*x\^2\+2\^2\*3\)/\(\(5\*x\^12\-2\*5\*x\^10\+\(2\^2\*5\+4\*2\*3\*5\+3\^2… | Problem multiple\_radicals\_0 \(vi\); inst a… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_1\-i | Trager | Sqrt\[x\^2\+1\]\*Sqrt\[1\-x\^4\] | Problem multiple\_radicals\_1 \(i\) | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_1\-ii | Trager | Sqrt\[a\^2\*x\^2\+b\^2\]\*Sqrt\[b\^4\-a\^4\*x\^4\] | Problem multiple\_radicals\_1 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| multrad\-multiple\_radicals\_1\-ii | PMT | Sqrt\[a\^2\*x\^2\+b\^2\]\*Sqrt\[b\^4\-a\^4\*x\^4\] | Problem multiple\_radicals\_1 \(ii\) | \{"failed", "tower construction failed", Sqrt\[b\^2 \+ a\^2 x\^2\] Sqrt\[b\^4 \- a\^4 x\^4\]\} |
| multrad\-multiple\_radicals\_1\-ii\-inst | Trager | Sqrt\[2\^2\*x\^2\+3\^2\]\*Sqrt\[3\^4\-2\^4\*x\^4\] | Problem multiple\_radicals\_1 \(ii\); inst a… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_1\-iii | Trager | 1/\(Sqrt\[x\^2\+1\]\*Sqrt\[1\-x\^4\]\) | Problem multiple\_radicals\_1 \(iii\) | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_1\-iv | Trager | 1/\(Sqrt\[a\^2\*x\^2\+b\^2\]\*Sqrt\[b\^4\-a\^4\*x\^4\]\) | Problem multiple\_radicals\_1 \(iv\) | \{"failed", "radicand is not a rational function over Q"\} |
| multrad\-multiple\_radicals\_1\-iv | PMT | 1/\(Sqrt\[a\^2\*x\^2\+b\^2\]\*Sqrt\[b\^4\-a\^4\*x\^4\]\) | Problem multiple\_radicals\_1 \(iv\) | \{"failed", "tower construction failed", 1/\(Sqrt\[b\^2 \+ a\^2 x\^2\] Sqrt\[b\^4 \- a\^4 x\^4\]\)\} |
| multrad\-multiple\_radicals\_1\-iv\-inst | Trager | 1/\(Sqrt\[2\^2\*x\^2\+3\^2\]\*Sqrt\[3\^4\-2\^4\*x\^4\]\) | Problem multiple\_radicals\_1 \(iv\); inst a… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_2\-i | Trager | \(x\^4\+2\)/\(\(x\^4\-2\)\*Sqrt\[x\^2\-1\]\*Sqrt\[x\^2\+2\]\) | Problem multiple\_radicals\_2 \(i\) | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_2\-ii | Trager | \(x\^4\-2\)\*Sqrt\[x\^2\+1\]\*Sqrt\[x\^2\+2\]/\(x\^2\*\(x\^4\+2\*x\^2\+2\)\) | Problem multiple\_radicals\_2 \(ii\) | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_2\-iii | Trager | \(x\^4\+2\*a\*x\^2\-a\*b\)\*Sqrt\[x\^2\+a\]/\(\(x\^4\-\(b\+1\)\*x\^2\-a\)\*\(x\^4\+\(1\-b\)\*x\^2\+a\)\*Sqrt\[x\^2\-b\]\) | Problem multiple\_radicals\_2 \(iii\) | \{"failed", "radicand is not a rational function over Q"\} |
| multrad\-multiple\_radicals\_2\-iii\-inst | Trager | \(x\^4\+2\*2\*x\^2\-2\*3\)\*Sqrt\[x\^2\+2\]/\(\(x\^4\-\(3\+1\)\*x\^2\-2\)\*\(x\^4\+\(1\-3\)\*x\^2\+2\)\*Sqrt\[x\^2\-3\]\) | Problem multiple\_radicals\_2 \(iii\); inst … | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_2\-iv | Trager | \(3\*x\^8\+\(a\-b\)\*x\^4\+a\*b\)\*Sqrt\[x\^4\+a\]\*Sqrt\[x\^4\-b\]/\(\(x\^8\+\(a\-b\)\*x\^4\+x\^2\-a\*b\)\*\(x\^8\+\(a\-b\)\*x\^4\-x\^2\-… | Problem multiple\_radicals\_2 \(iv\) | \{"failed", "radicand is not a rational function over Q"\} |
| multrad\-multiple\_radicals\_2\-iv | PMT | \(3\*x\^8\+\(a\-b\)\*x\^4\+a\*b\)\*Sqrt\[x\^4\+a\]\*Sqrt\[x\^4\-b\]/\(\(x\^8\+\(a\-b\)\*x\^4\+x\^2\-a\*b\)\*\(x\^8\+\(a\-b\)\*x\^4\-x\^2\-… | Problem multiple\_radicals\_2 \(iv\) | \{"failed", "tower construction failed", \(Sqrt\[a \+ x\^4\] Sqrt\[\-b \+ x\^4\] \(a b \+ \(a \- b\) x\^4 \+ 3 x\^8\)\)/\(\(\-a b \- x\^… |
| multrad\-multiple\_radicals\_2\-iv\-inst | Trager | \(3\*x\^8\+\(2\-3\)\*x\^4\+2\*3\)\*Sqrt\[x\^4\+2\]\*Sqrt\[x\^4\-3\]/\(\(x\^8\+\(2\-3\)\*x\^4\+x\^2\-2\*3\)\*\(x\^8\+\(2\-3\)\*x\^4\-x\^2\-… | Problem multiple\_radicals\_2 \(iv\); inst a… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_2\-iv\-inst | PMT | \(3\*x\^8\+\(2\-3\)\*x\^4\+2\*3\)\*Sqrt\[x\^4\+2\]\*Sqrt\[x\^4\-3\]/\(\(x\^8\+\(2\-3\)\*x\^4\+x\^2\-2\*3\)\*\(x\^8\+\(2\-3\)\*x\^4\-x\^2\-… | Problem multiple\_radicals\_2 \(iv\); inst a… | \{"failed", "tower construction failed", \(Sqrt\[\-3 \+ x\^4\] Sqrt\[2 \+ x\^4\] \(6 \- x\^4 \+ 3 x\^8\)\)/\(\(\-6 \- x\^2 \- x\^4 \+ x\^… |
| multrad\-multiple\_radicals\_2\-v\-inst | Trager | \(2\*x\^5\-3\*3\*x\^4\+4\*5\*x\^3\-2\*2\*x\+2\*3\)\*Sqrt\[x\^4\+2\]/\(\(x\^4\-x\^2\+3\*x\+2\-5\)\*\(x\^2\-3\*x\+5\)\^\(3/2\)\) | Problem multiple\_radicals\_2 \(v\); inst a=… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_2\-v\-inst | PMT | \(2\*x\^5\-3\*3\*x\^4\+4\*5\*x\^3\-2\*2\*x\+2\*3\)\*Sqrt\[x\^4\+2\]/\(\(x\^4\-x\^2\+3\*x\+2\-5\)\*\(x\^2\-3\*x\+5\)\^\(3/2\)\) | Problem multiple\_radicals\_2 \(v\); inst a=… | \{"failed", "tower construction failed", \(Sqrt\[2 \+ x\^4\] \(6 \- 4 x \+ 20 x\^3 \- 9 x\^4 \+ 2 x\^5\)\)/\(\(5 \- 3 x \+ x\^2\)\^\(3… |
| multrad\-multiple\_radicals\_3\-i\-inst | Trager | \(x\^2\-2\*2\*x\+2\*3\+2\*5\-3\*5\)\*Sqrt\[x\-3\]\*Sqrt\[x\-5\]/\(\(x\^2\-\(3\+5\-1\)\*x\+3\*5\-2\)\*\(x\-2\)\^\(3/2\)\) | Problem multiple\_radicals\_3 \(i\); inst a=… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_3\-ii | Trager | \(x\^2\-2\*a\*x\+a\*b\+a\*c\-b\*c\)/\(\(x\^2\-\(b\+c\+1\)\*x\+a\+b\*c\)\*Sqrt\[x\-a\]\*Sqrt\[x\-b\]\*Sqrt\[x\-c\]\) | Problem multiple\_radicals\_3 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| multrad\-multiple\_radicals\_3\-ii\-inst | Trager | \(x\^2\-2\*2\*x\+2\*3\+2\*5\-3\*5\)/\(\(x\^2\-\(3\+5\+1\)\*x\+2\+3\*5\)\*Sqrt\[x\-2\]\*Sqrt\[x\-3\]\*Sqrt\[x\-5\]\) | Problem multiple\_radicals\_3 \(ii\); inst a… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_3\-iii | Trager | \(x\^2\-2\*b\*x\+a\*b\-a\*c\+b\*c\)\*Sqrt\[x\-b\]/\(\(x\^2\-\(a\+c\-1\)\*x\-b\+a\*c\)\*\(x\^2\-\(a\+c\+1\)\*x\+b\+a\*c\)\*Sqrt\[x\-a\]\*S… | Problem multiple\_radicals\_3 \(iii\) | \{"failed", "radicand is not a rational function over Q"\} |
| multrad\-multiple\_radicals\_3\-iii\-inst | Trager | \(x\^2\-2\*3\*x\+2\*3\-2\*5\+3\*5\)\*Sqrt\[x\-3\]/\(\(x\^2\-\(2\+5\-1\)\*x\-3\+2\*5\)\*\(x\^2\-\(2\+5\+1\)\*x\+3\+2\*5\)\*Sqrt\[x\-2\]\*S… | Problem multiple\_radicals\_3 \(iii\); inst … | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_4\-i | Trager | Sqrt\[1\-x\]\*Sqrt\[x\+1\]\*Sqrt\[1\-x\^4\] | Problem multiple\_radicals\_4 \(i\) | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_4\-ii | Trager | Sqrt\[b\-a\*x\]\*Sqrt\[a\*x\+b\]\*Sqrt\[b\^4\-a\^4\*x\^4\] | Problem multiple\_radicals\_4 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| multrad\-multiple\_radicals\_4\-ii | PMT | Sqrt\[b\-a\*x\]\*Sqrt\[a\*x\+b\]\*Sqrt\[b\^4\-a\^4\*x\^4\] | Problem multiple\_radicals\_4 \(ii\) | \{"failed", "tower construction failed", Sqrt\[b \- a x\] Sqrt\[b \+ a x\] Sqrt\[b\^4 \- a\^4 x\^4\]\} |
| multrad\-multiple\_radicals\_4\-ii\-inst | Trager | Sqrt\[3\-2\*x\]\*Sqrt\[2\*x\+3\]\*Sqrt\[3\^4\-2\^4\*x\^4\] | Problem multiple\_radicals\_4 \(ii\); inst a… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_4\-ii\-inst | PMT | Sqrt\[3\-2\*x\]\*Sqrt\[2\*x\+3\]\*Sqrt\[3\^4\-2\^4\*x\^4\] | Problem multiple\_radicals\_4 \(ii\); inst a… | \{"failed", "tower construction failed", Sqrt\[3 \- 2 x\] Sqrt\[3 \+ 2 x\] Sqrt\[81 \- 16 x\^4\]\} |
| multrad\-multiple\_radicals\_4\-iii | Trager | Sqrt\[x\+1\]\*Sqrt\[1\-x\^4\]/Sqrt\[1\-x\] | Problem multiple\_radicals\_4 \(iii\) | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_4\-iv | Trager | Sqrt\[a\*x\+b\]\*Sqrt\[b\^4\-a\^4\*x\^4\]/Sqrt\[b\-a\*x\] | Problem multiple\_radicals\_4 \(iv\) | \{"failed", "radicand is not a rational function over Q"\} |
| multrad\-multiple\_radicals\_4\-iv | PMT | Sqrt\[a\*x\+b\]\*Sqrt\[b\^4\-a\^4\*x\^4\]/Sqrt\[b\-a\*x\] | Problem multiple\_radicals\_4 \(iv\) | \{"failed", "tower construction failed", \(Sqrt\[b \+ a x\] Sqrt\[b\^4 \- a\^4 x\^4\]\)/Sqrt\[b \- a x\]\} |
| multrad\-multiple\_radicals\_4\-iv\-inst | Trager | Sqrt\[2\*x\+3\]\*Sqrt\[3\^4\-2\^4\*x\^4\]/Sqrt\[3\-2\*x\] | Problem multiple\_radicals\_4 \(iv\); inst a… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_4\-iv\-inst | PMT | Sqrt\[2\*x\+3\]\*Sqrt\[3\^4\-2\^4\*x\^4\]/Sqrt\[3\-2\*x\] | Problem multiple\_radicals\_4 \(iv\); inst a… | \{"failed", "tower construction failed", \(Sqrt\[3 \+ 2 x\] Sqrt\[81 \- 16 x\^4\]\)/Sqrt\[3 \- 2 x\]\} |
| multrad\-multiple\_radicals\_4\-v | Trager | Sqrt\[x\+1\]/\(Sqrt\[1\-x\]\*Sqrt\[1\-x\^4\]\) | Problem multiple\_radicals\_4 \(v\) | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_4\-v | PMT | Sqrt\[x\+1\]/\(Sqrt\[1\-x\]\*Sqrt\[1\-x\^4\]\) | Problem multiple\_radicals\_4 \(v\) | \{"needs torsion realisation \(milestone iii\)", \{\{\-1 \+ \(\-1 \+ Sqrt\[1 \- x\]\)\^2/\(\-1 \+ Sqrt\[1 \+ x\]\)\^2 \+ \(2 \(\-1 \+ Sqrt… |
| multrad\-multiple\_radicals\_4\-vi | Trager | Sqrt\[a\*x\+b\]/\(Sqrt\[b\-a\*x\]\*Sqrt\[b\^4\-a\^4\*x\^4\]\) | Problem multiple\_radicals\_4 \(vi\) | \{"failed", "radicand is not a rational function over Q"\} |
| multrad\-multiple\_radicals\_4\-vi | PMT | Sqrt\[a\*x\+b\]/\(Sqrt\[b\-a\*x\]\*Sqrt\[b\^4\-a\^4\*x\^4\]\) | Problem multiple\_radicals\_4 \(vi\) | \{"failed", "tower construction failed", Sqrt\[b \+ a x\]/\(Sqrt\[b \- a x\] Sqrt\[b\^4 \- a\^4 x\^4\]\)\} |
| multrad\-multiple\_radicals\_4\-vi\-inst | Trager | Sqrt\[2\*x\+3\]/\(Sqrt\[3\-2\*x\]\*Sqrt\[3\^4\-2\^4\*x\^4\]\) | Problem multiple\_radicals\_4 \(vi\); inst a… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_4\-vi\-inst | PMT | Sqrt\[2\*x\+3\]/\(Sqrt\[3\-2\*x\]\*Sqrt\[3\^4\-2\^4\*x\^4\]\) | Problem multiple\_radicals\_4 \(vi\); inst a… | \{"failed", "tower construction failed", Sqrt\[3 \+ 2 x\]/\(Sqrt\[3 \- 2 x\] Sqrt\[81 \- 16 x\^4\]\)\} |
| multrad\-multiple\_radicals\_4\-vii | Trager | 1/\(Sqrt\[1\-x\]\*Sqrt\[x\+1\]\*Sqrt\[1\-x\^4\]\) | Problem multiple\_radicals\_4 \(vii\) | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_4\-vii | PMT | 1/\(Sqrt\[1\-x\]\*Sqrt\[x\+1\]\*Sqrt\[1\-x\^4\]\) | Problem multiple\_radicals\_4 \(vii\) | \{"needs torsion realisation \(milestone iii\)", \{\{\-1 \+ \(\-1 \+ Sqrt\[1 \- x\]\)\^2/\(\-1 \+ Sqrt\[1 \+ x\]\)\^2 \- \(2 \(\-1 \+ Sqrt… |
| multrad\-multiple\_radicals\_4\-viii | Trager | 1/\(Sqrt\[b\-a\*x\]\*Sqrt\[a\*x\+b\]\*Sqrt\[b\^4\-a\^4\*x\^4\]\) | Problem multiple\_radicals\_4 \(viii\) | \{"failed", "radicand is not a rational function over Q"\} |
| multrad\-multiple\_radicals\_4\-viii | PMT | 1/\(Sqrt\[b\-a\*x\]\*Sqrt\[a\*x\+b\]\*Sqrt\[b\^4\-a\^4\*x\^4\]\) | Problem multiple\_radicals\_4 \(viii\) | \{"failed", "tower construction failed", 1/\(Sqrt\[b \- a x\] Sqrt\[b \+ a x\] Sqrt\[b\^4 \- a\^4 x\^4\]\)\} |
| multrad\-multiple\_radicals\_4\-viii\-inst | Trager | 1/\(Sqrt\[3\-2\*x\]\*Sqrt\[2\*x\+3\]\*Sqrt\[3\^4\-2\^4\*x\^4\]\) | Problem multiple\_radicals\_4 \(viii\); inst… | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| multrad\-multiple\_radicals\_4\-viii\-inst | PMT | 1/\(Sqrt\[3\-2\*x\]\*Sqrt\[2\*x\+3\]\*Sqrt\[3\^4\-2\^4\*x\^4\]\) | Problem multiple\_radicals\_4 \(viii\); inst… | \{"failed", "tower construction failed", 1/\(Sqrt\[3 \- 2 x\] Sqrt\[3 \+ 2 x\] Sqrt\[81 \- 16 x\^4\]\)\} |
| multrad\-multiple\_radicals\_7\-ii | Trager | ArcTan\[Sqrt\[x\+1\]\-Sqrt\[x\]\] | Problem multiple\_radicals\_7 \(ii\), Charlw… | \{"failed", "transcendental function in the integrand"\} |
| multrad\-multiple\_radicals\_7\-iv | Trager | ArcTan\[d\*Sqrt\[a\*x\+b\]\-d\*Sqrt\[a\*x\+c\]\] | Problem multiple\_radicals\_7 \(iv\) | \{"failed", "transcendental function in the integrand"\} |
| multrad\-multiple\_radicals\_7\-iv\-inst | Trager | ArcTan\[7\*Sqrt\[2\*x\+3\]\-7\*Sqrt\[2\*x\+5\]\] | Problem multiple\_radicals\_7 \(iv\); inst a… | \{"failed", "transcendental function in the integrand"\} |
| euler\-prob\_euler\_outrageous\_1\-ii | Trager | \(a\^2\*x\^2\-2\*a\*b\*x\+b\^2\)\*\(a\^2\*x\^2\+2\*a\*b\*x\+b\^2\)/\(\(a\^2\*x\^2\+b\^2\)\*\(a\^4\*x\^4\+6\*a\^2\*b\^2\*x\^2\+b\^4\)\^\(3/… | Problem prob\_euler\_outrageous\_1 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| euler\-prob\_euler\_outrageous\_2\-p | PMT | 1/\(\(3\+x\^2\)\*\(1\+3\*x\^2\)\^\(1/3\)\) | Problem prob\_euler\_outrageous\_2 \(\+ sign\)… | \{"failed", "residue class lost: no places over the prime for a root of the residue polynomial", 3 \+ x\^2, \(\-1/4… |
| euler\-prob\_euler\_outrageous\_3\-ii | Trager | \(a\*x\^2\-x\-2\*b\)/\(\(x\^4\-a\*x\^2\+b\)\*\(a\*x\^2\-b\)\^\(1/4\)\) | Problem prob\_euler\_outrageous\_3 \(ii\) | \{"failed", "radicand is not a rational function over Q"\} |
| euler\-prob\_euler\_outrageous\_3\-ii | PMT | \(a\*x\^2\-x\-2\*b\)/\(\(x\^4\-a\*x\^2\+b\)\*\(a\*x\^2\-b\)\^\(1/4\)\) | Problem prob\_euler\_outrageous\_3 \(ii\) | \{"failed", "residue class lost: no places over the prime for a root of the residue polynomial", \-b \+ a x\^2 \- x… |
| euler\-prob\_euler\_outrageous\_3\-ii\-inst | PMT | \(2\*x\^2\-x\-2\*3\)/\(\(x\^4\-2\*x\^2\+3\)\*\(2\*x\^2\-3\)\^\(1/4\)\) | Problem prob\_euler\_outrageous\_3 \(ii\); in… | \{"failed", "residue class lost: no places over the prime for a root of the residue polynomial", 3 \- 2 x\^2 \+ x\^… |
| euler\-prob\_euler\_outrageous\_4\-iii | Trager | Sqrt\[a\*x\^4\+c\*x\^2\+b\]/\(b\-a\*x\^4\) | Problem prob\_euler\_outrageous\_4 \(iii\) | \{"failed", "radicand is not a rational function over Q"\} |
| euler\-prob\_euler\_outrageous\_4\-iii\-inst | PMT | Sqrt\[2\*x\^4\+5\*x\^2\+3\]/\(3\-2\*x\^4\) | Problem prob\_euler\_outrageous\_4 \(iii\); i… | \{"needs torsion realisation \(milestone iii\)", \{\{\-3 \+ 2 x\^4, \{Root\[1 \- 960 \#1\^2 \+ 9216 \#1\^4 \&, 2\], Root\[1 \- 960… |
| euler\-prob\_euler\_outrageous\_4\-iv | Trager | x\^2/\(\(b\-a\*x\^4\)\*Sqrt\[a\*x\^4\+b\]\) | Problem prob\_euler\_outrageous\_4 \(iv\) | \{"failed", "radicand is not a rational function over Q"\} |
| euler\-prob\_euler\_outrageous\_5 | Trager | \(a\*x\^2\-b\)/\(\(a\*x\^2\+b\)\*Sqrt\[a\^2\*x\^4\+b\^2\]\) | Problem prob\_euler\_outrageous\_5 | \{"failed", "radicand is not a rational function over Q"\} |
| euler\-prob\_euler\_outrageous\_6\-i | Trager | 1/\(\(x\^4\+1\)\*Sqrt\[Sqrt\[x\^4\+1\]\-x\^2\]\) | Problem prob\_euler\_outrageous\_6 \(i\) | \{"failed", "nested radicals"\} |
| euler\-prob\_euler\_outrageous\_6\-i | PMT | 1/\(\(x\^4\+1\)\*Sqrt\[Sqrt\[x\^4\+1\]\-x\^2\]\) | Problem prob\_euler\_outrageous\_6 \(i\) | \{"failed", "tower construction failed", 1/\(\(1 \+ x\^4\) Sqrt\[\-x\^2 \+ Sqrt\[1 \+ x\^4\]\]\)\} |
| euler\-prob\_euler\_outrageous\_6\-ii | Trager | 1/\(\(a\^2\*x\^4\+b\)\*Sqrt\[Sqrt\[a\^2\*x\^4\+b\]\-a\*x\^2\]\) | Problem prob\_euler\_outrageous\_6 \(ii\) | \{"failed", "nested radicals"\} |
| euler\-prob\_euler\_outrageous\_6\-ii | PMT | 1/\(\(a\^2\*x\^4\+b\)\*Sqrt\[Sqrt\[a\^2\*x\^4\+b\]\-a\*x\^2\]\) | Problem prob\_euler\_outrageous\_6 \(ii\) | \{"failed", "tower construction failed", 1/\(\(b \+ a\^2 x\^4\) Sqrt\[\-a x\^2 \+ Sqrt\[b \+ a\^2 x\^4\]\]\)\} |
| euler\-prob\_euler\_outrageous\_6\-ii\-inst | Trager | 1/\(\(2\^2\*x\^4\+3\)\*Sqrt\[Sqrt\[2\^2\*x\^4\+3\]\-2\*x\^2\]\) | Problem prob\_euler\_outrageous\_6 \(ii\); in… | \{"failed", "nested radicals"\} |
| euler\-prob\_euler\_outrageous\_6\-ii\-inst | PMT | 1/\(\(2\^2\*x\^4\+3\)\*Sqrt\[Sqrt\[2\^2\*x\^4\+3\]\-2\*x\^2\]\) | Problem prob\_euler\_outrageous\_6 \(ii\); in… | \{"failed", "tower construction failed", 1/\(\(3 \+ 4 x\^4\) Sqrt\[\-2 x\^2 \+ Sqrt\[3 \+ 4 x\^4\]\]\)\} |
| euler\-prob\_euler\_outrageous\_6\-iii\-n1 | Trager | 1/\(\(x\^2\+1\)\*Sqrt\[\(x\^2\+1\)\-x\^2\]\) | Problem prob\_euler\_outrageous\_6 \(iii\), n… | \{"failed", "no radical in the integrand \(a rational function\)"\} |
| euler\-prob\_euler\_outrageous\_6\-iii\-n2 | Trager | 1/\(\(x\^4\+1\)\*Sqrt\[\(x\^4\+1\)\^\(1/2\)\-x\^2\]\) | Problem prob\_euler\_outrageous\_6 \(iii\), n… | \{"failed", "nested radicals"\} |
| euler\-prob\_euler\_outrageous\_6\-iii\-n2 | PMT | 1/\(\(x\^4\+1\)\*Sqrt\[\(x\^4\+1\)\^\(1/2\)\-x\^2\]\) | Problem prob\_euler\_outrageous\_6 \(iii\), n… | \{"failed", "tower construction failed", 1/\(\(1 \+ x\^4\) Sqrt\[\-x\^2 \+ Sqrt\[1 \+ x\^4\]\]\)\} |
| euler\-prob\_euler\_outrageous\_6\-iii\-n3 | Trager | 1/\(\(x\^6\+1\)\*Sqrt\[\(x\^6\+1\)\^\(1/3\)\-x\^2\]\) | Problem prob\_euler\_outrageous\_6 \(iii\), n… | \{"failed", "nested radicals"\} |
| euler\-prob\_euler\_outrageous\_6\-iii\-n3 | PMT | 1/\(\(x\^6\+1\)\*Sqrt\[\(x\^6\+1\)\^\(1/3\)\-x\^2\]\) | Problem prob\_euler\_outrageous\_6 \(iii\), n… | \{"failed", "tower construction failed", 1/\(\(1 \+ x\^6\) Sqrt\[\-x\^2 \+ \(1 \+ x\^6\)\^\(1/3\)\]\)\} |
| timofeev\-timofeev\_1\-i | PMT | 1/\(\(x\+1\)\^2\*\(x\-1\)\^4\)\^\(1/3\) | Problem timofeev\_1 \(i\) | \{"failed", "verification failed"\} |
| timofeev\-timofeev\_3\-i | Trager | Sqrt\[x\-5\]\*Sqrt\[x\+3\]/\(\(x\-1\)\*\(x\^2\-25\)\) | Problem timofeev\_3 \(i\) | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| timofeev\-timofeev\_3\-ii | Trager | x\^2\*Sqrt\[x\+1\]\*\(1\-x\^2\)\^\(1/4\)/\(Sqrt\[1\-x\]\*\(Sqrt\[1\-x\]\-Sqrt\[x\+1\]\)\) | Problem timofeev\_3 \(ii\) | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| timofeev\-timofeev\_3\-ii | PMT | x\^2\*Sqrt\[x\+1\]\*\(1\-x\^2\)\^\(1/4\)/\(Sqrt\[1\-x\]\*\(Sqrt\[1\-x\]\-Sqrt\[x\+1\]\)\) | Problem timofeev\_3 \(ii\) | \{"failed", "partially ramified branch prime", \-1 \+ \(\-1 \+ Sqrt\[1 \- x\]\)\^2/\(\-1 \+ Sqrt\[1 \+ x\]\)\^2 \- \(2 \(\-1 \+ Sqrt\[1… |
| timofeev\-timofeev\_3\-iii | Trager | x\*Sqrt\[1\-x\]\*\(x\+1\)\^\(2/3\)/\(\(1\-x\)\^\(2/3\)\*Sqrt\[x\+1\]\-\(1\-x\)\^\(5/6\)\*\(x\+1\)\^\(1/3\)\) | Problem timofeev\_3 \(iii\) | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| timofeev\-timofeev\_3\-iii | PMT | x\*Sqrt\[1\-x\]\*\(x\+1\)\^\(2/3\)/\(\(1\-x\)\^\(2/3\)\*Sqrt\[x\+1\]\-\(1\-x\)\^\(5/6\)\*\(x\+1\)\^\(1/3\)\) | Problem timofeev\_3 \(iii\) | \{"failed", "no solution within bounds", \{12\}\} |
| timofeev\-timofeev\_3\-v | Trager | \(\(x\^2\+2\)/x\^2\)\^\(7/9\)/\(x\^2\+2\)\^\(3/2\) | Problem timofeev\_3 \(v\) | \{"failed", "more than one radicand \(a radical compositum\)"\} |
| potpourri\-p1 | PMT | \(x\^2\+1\)\*\(x\^2\-x\-1\)/\(\(x\^2\-1\)\^2\*Sqrt\[x\^4\-x\^3\-x\^2\+x\+1\]\) | Potpourri problem 1 | \{"needs torsion realisation \(milestone iii\)", \{\{\-1 \+ x, \{1/2, \-1/2\}\}, \{1 \+ x, \{\-1/2, 1/2\}\}\}\} |
| potpourri\-p2 | PMT | 1/\(\(x\-1\)\*Sqrt\[3\*x\^4\-12\*x\^3\+18\*x\^2\-12\*x\-1\]\) | Potpourri problem 2 | \{"failed", "no solution within bounds", \{0\}\} |
| potpourri\-p6\-vii | PMT | 1/\(x\^8\+2\*x\^4\+1\)\^\(1/8\) | Potpourri problem 6 \(vii\) | \{"failed", "partially ramified branch prime", 1 \+ x\^4\} |
| potpourri\-p6\-viii | PMT | \(x\^3\-1\)\*Sqrt\[x\^6\-1\]/\(x\^4\*\(x\^3\+1\)\) | Potpourri problem 6 \(viii\) | \{"needs torsion realisation \(milestone iii\)", \{\{x, \{2\*I, \-2\*I\}\}\}\} |
| potpourri\-p6\-ix | PMT | \(2\*x\^4\-1\)/\(\(x\^8\-1\)\*\(x\^4\-1\)\^\(1/4\)\) | Potpourri problem 6 \(ix\) | \{"failed", "no solution within bounds", \{4\}\} |
| potpourri\-p8\-iv | Trager | x/\(\(x\^2\-1\)\*Sqrt\[\-x\+Sqrt\[x\^2\+1\]\]\) | Potpourri problem 8 \(nested radicals\) \(i… | \{"failed", "nested radicals"\} |
| potpourri\-p8\-vi | Trager | Sqrt\[1\+Sqrt\[x\+Sqrt\[x\^2\+1\]\]\]/Sqrt\[x\^2\+1\] | Potpourri problem 8 \(nested radicals\) \(v… | \{"failed", "nested radicals"\} |
| potpourri\-p8\-viii | Trager | \(1\+\(x\+Sqrt\[x\^2\-1\]\)\^\(1/4\)\)\^\(1/6\)/Sqrt\[x\^2\-1\] | Potpourri problem 8 \(nested radicals\) \(v… | \{"failed", "nested radicals"\} |
| potpourri\-p9\-i | Trager | \(a\+Sqrt\[x\]\+Sqrt\[x\+2\*a\*Sqrt\[x\]\]\)\^\(1/3\) | Potpourri problem 9 \(parametric nested r… | \{"failed", "nested radicals"\} |
| potpourri\-p9\-i\-inst | Trager | \(2\+Sqrt\[x\]\+Sqrt\[x\+2\*2\*Sqrt\[x\]\]\)\^\(1/3\) | Potpourri problem 9 \(parametric nested r… | \{"failed", "nested radicals"\} |
| potpourri\-p9\-ii | Trager | Sqrt\[c\*x\^2\-x\*Sqrt\[a\*x\^2\-b\*x\]\]/\(x\^4\*Sqrt\[a\*x\^2\-b\*x\]\) | Potpourri problem 9 \(parametric nested r… | \{"failed", "nested radicals"\} |
| potpourri\-p9\-ii | PMT | Sqrt\[c\*x\^2\-x\*Sqrt\[a\*x\^2\-b\*x\]\]/\(x\^4\*Sqrt\[a\*x\^2\-b\*x\]\) | Potpourri problem 9 \(parametric nested r… | \{"failed", "tower construction failed", Sqrt\[c x\^2 \- x Sqrt\[\-b x \+ a x\^2\]\]/\(x\^4 Sqrt\[\-b x \+ a x\^2\]\)\} |
| potpourri\-p9\-ii\-inst | Trager | Sqrt\[5\*x\^2\-x\*Sqrt\[2\*x\^2\-3\*x\]\]/\(x\^4\*Sqrt\[2\*x\^2\-3\*x\]\) | Potpourri problem 9 \(parametric nested r… | \{"failed", "nested radicals"\} |
| potpourri\-p9\-ii\-inst | PMT | Sqrt\[5\*x\^2\-x\*Sqrt\[2\*x\^2\-3\*x\]\]/\(x\^4\*Sqrt\[2\*x\^2\-3\*x\]\) | Potpourri problem 9 \(parametric nested r… | \{"failed", "tower construction failed", Sqrt\[5 x\^2 \- x Sqrt\[\-3 x \+ 2 x\^2\]\]/\(x\^4 Sqrt\[\-3 x \+ 2 x\^2\]\)\} |
| danger\-d3\-i | PMT | 1/\(x\^4\-4\*x\^2\+2\)\^\(1/4\) | Danger zone problem 3; todo says not tha… | \{"failed", "no solution within bounds", \{0\}\} |
| danger\-d3\-ii | PMT | 1/\(x\^4\+4\*x\^2\+2\)\^\(1/4\) | Danger zone problem 3; todo says not tha… | \{"failed", "no solution within bounds", \{0\}\} |
| danger\-d3\-iii | PMT | 1/\(x\^4\-8\*x\^2\+8\)\^\(1/4\) | Danger zone problem 3; todo says not tha… | \{"failed", "no solution within bounds", \{0\}\} |
| danger\-d3\-iv | PMT | 1/\(x\^4\+8\*x\^2\+8\)\^\(1/4\) | Danger zone problem 3; todo says not tha… | \{"failed", "no solution within bounds", \{0\}\} |
| danger\-d4\-i | PMT | 1/\(x\^4\-4\*x\^3\+2\*x\^2\)\^\(1/4\) | Danger zone problem 4; todo says not tha… | \{"failed", "partially ramified branch prime", x\} |
| danger\-d4\-ii | PMT | 1/\(x\^4\+4\*x\^3\+2\*x\^2\)\^\(1/4\) | Danger zone problem 4; todo says not tha… | \{"failed", "partially ramified branch prime", x\} |
| danger\-d4\-iii | PMT | 1/\(x\^4\-8\*x\^3\+8\*x\^2\)\^\(1/4\) | Danger zone problem 4; todo says not tha… | \{"failed", "partially ramified branch prime", x\} |
| danger\-d4\-iv | PMT | 1/\(x\^4\+8\*x\^3\+8\*x\^2\)\^\(1/4\) | Danger zone problem 4; todo says not tha… | \{"failed", "partially ramified branch prime", x\} |
| intro\-can\_cannot\-1 | Trager | Sqrt\[Tan\[x\]\] | Introduction, Integrals that Can and Can… | \{"failed", "transcendental function in the integrand"\} |
| intro\-can\_cannot\-2 | Trager | ArcTan\[Sqrt\[x\]\] | Introduction, Integrals that Can and Can… | \{"failed", "transcendental function in the integrand"\} |
| intro\-can\_cannot\-3 | Trager | Tan\[Sqrt\[x\]\]/Sqrt\[x\] | Introduction, Integrals that Can and Can… | \{"failed", "transcendental function in the integrand"\} |
| intro\-can\_cannot\-4 | Trager | ArcTan\[x\]/Sqrt\[x\] | Introduction, Integrals that Can and Can… | \{"failed", "transcendental function in the integrand"\} |
| intro\-validity\-2 | Trager | \(x\^8\-1\)/\(\(x\^8\-a\*x\^4\+1\)\*Sqrt\[x\^4\-1\]\) | Introduction, Validity of Solutions, par… | \{"failed", "integrand is not a rational function of x and the radical over Q"\} |
| intro\-validity\-2 | PMT | \(x\^8\-1\)/\(\(x\^8\-a\*x\^4\+1\)\*Sqrt\[x\^4\-1\]\) | Introduction, Validity of Solutions, par… | \{"failed", "residue class lost: no places over the prime for a root of the residue polynomial", \-1 \+ a x\^4 \- x… |
| intro\-log\-1 | Trager | 1/x | Introduction, integral of 1/x section | \{"failed", "no radical in the integrand \(a rational function\)"\} |
| intro\-multiple\_approaches\-1 | Trager | 1/\(\(x\-1\)\*Sqrt\[x\^2\+b\*x\+c\]\) | Introduction, Multiple Approaches to Int… | \{"failed", "radicand is not a rational function over Q"\} |

### FriCAS answers that fail their own derivative check (1)

- `potpourri-p7-iii`: `1/(x^6+2*x^5-x^4-4*x^3-x^2+2*x+1)^(1/6)`

<!-- END GENERATED -->
