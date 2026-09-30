# `Integrate`ParallelMixedSpecial`: current state, and the route to 305/312

The special-function stage of S. Blake's parallel Risch–Norman integrator landed
in **v0.239** (`src/internal/ParallelMixedSpecial.m`, 1868 lines, on top of the
Part II delta in `src/internal/mixed/ParallelMixed.m`). This document records
exactly where it stands against the three reference ports and what each remaining
case needs, so the work can be picked up without re-deriving any of it.

Measurement: `<research>/special/stress_mathilda.py`, 312 cases, `TIMEOUT=120`,
`WORKERS=6`, verdicts scored by `stress_judge.py` — the same rule all four ports
are scored by, not a re-implementation. Re-run with
`python3 stress_compare.py --ports=python,wl,maxima,mathilda`.

---

## 1. Where it stands

| port | PASS | HONEST | WEAK | FAIL |
|---|---|---|---|---|
| Python (SymPy, the paper's reference) | 305 | 7 | 0 | 0 |
| Mathematica 14 | 305 | 7 | 0 | 0 |
| Maxima | 305 | 7 | 0 | 0 |
| **Mathilda v0.239** | **247** | **51** | **11** | **3** |

**No false certificate anywhere.** That is the bar this corpus exists to police —
a `{"not elementary", …}` for an integrand that *is* elementary is the one
failure mode worse than declining — and Mathilda clears it.

Seven of the seventeen groups are already at parity: G not-in-class 14/14,
H partial 10/10, Q random partial 20/20, R partial (elementary) 15/15,
R partial (radicals) 5/5, S random (elementary) 15/15, and E radical (the same
2 HONEST as every port). The elliptic group F is 13/15.

The answers themselves are a faithful port. Where Mathilda withholds, the answer
it computed is **character-for-character Mathematica's**:

| case | Mathilda's computed answer | reference |
|---|---|---|
| `exp(-x^2)` | `-Sqrt[Pi] Erfc[x]/2` | identical |
| `exp(-x^3)` | `-Gamma[1/3, x^3]/3` | identical |
| `log(x) sin(x)` | `CosIntegral[x] - Log[x]/(1 + Tan[x/2]^2) + (Log[x] Tan[x/2]^2)/(1 + Tan[x/2]^2)` | identical |

So the gap is **not** in the kernels, the ansatz, the bounds, or the surface
reconstruction. It is almost entirely in one place: the non-elementarity
certificate that strict mode requires before it will hand back a special answer.

### The 58 differing cases, bucketed by root cause

| # | cause | cases | ids |
|---|---|---|---|
| **A** | certificate withheld — the answer is found and correct, strict mode suppresses it | **39** | 55 56 65–83 94 157 159 161 163–168 170 271–273 278 279 308 312 |
| **B** | partial remainder — the partial mode left exactly the sub-integrands of (A) in the remainder | **11** | 226 230 231 233 235 250 252 258 264 266 270 |
| **C** | no solution with the special kernels — genuine coverage | **5** | 87 89 90 172 175 |
| **D** | timeout at 120 s on the mod-p non-torsion certificate | **2** | 109 110 |
| **E** | unverified answer — branch fidelity | **1** | 34 |

**A + B = 50 of 58 share one cause.** Closing it is the whole of the work between
247 and ~297; C, D and E are three separate, much smaller jobs.

(Group counts for A: C erf/Gamma 19, P families 10, M Cherry 5, B Ei/li 2,
K structure 2, D trigonometric 1. For B: S random (special) 6, R partial
(special) 5.)

---

## 2. Work item 1 — the certificate (50 cases)

### What is actually happening

`ParallelIntegrateSpecial` (Algorithm S7, `src/internal/ParallelMixedSpecial.m`)
finds the special answer, then asks for a certificate of non-elementarity before
returning it — Theorem 8.6, and the reason the stage is sound. With
`$StrictSP = True` (the default, and what all four ports ran) no certificate means

```
{"failed", "special answer found, but the integrand is not certified non-elementary", answer}
```

The certificate comes from `CertifyNonelementary[f, T]`, which asks Part II
(`Part2[f, T, ext]`, for `ext` in `{False, True}`, over the tower and over its
pruned tower) for a `{"not elementary", …}`. For every case in bucket A, Part II
answers `{"failed", "no solution within bounds", …}` instead.

### Two findings, both measured, that redirect the obvious fix

**Finding 1 — the `q =!= None` guard is expired but is not the blocker.**

`src/internal/mixed/ParallelMixed.m` guards the holomorphic-remainder
certificate with `q =!= None`; the research reference has no such guard. It was
added at v0.161 as "the T2/T10 false-certificate protection" because Mathilda
emitted **false** certificates on curve-free towers for

```
T2   1/(x (Log[x]^2 + 1))    -> ArcTan[Log[x]]     (elementary)
T10  1/(x^4 - 1)             -> elementary
```

`ArcTan` is a logarithm over `Q(i)`, and the defect that made Mathilda miss it
was `FreeQ[e, Complex]` returning True on a Gaussian atom (A19), so `I` never
entered the field and the logand set was incomplete — an incomplete logand set
makes the system inconsistent and the certificate fire falsely. **That defect was
fixed at v0.175** (the AlgAtoms fix). Measured: with the guard lifted behind a
temporary flag, **T2 and T10 both solve correctly** (`ArcTan[Log[x]]` and the
rational answer). So the guard's original justification has expired.

But lifting it is **not sufficient**: with the guard lifted, `Exp[-x^2]`,
`Exp[x^2]` and `Sqrt[x] Exp[x]` still answer `{"failed", "no solution within
bounds", {2, 1}}`. Something downstream of the guard blocks the certificate.

**Finding 2 — `VerifiedResidueFree` is order-dependent, and conflates "the probe
did not finish" with "the residual has residues".**

The certificate's last conjunct is `VerifiedResidueFree[rem, T]`
(`src/internal/mixed/ParallelMixed.m`):

```mathematica
VerifiedResidueFree[rem_, T_] := Module[{key = Hash[{rem, T}], ...},
  Catch[iPIM[rem, T, "Bounds" -> ConstantArray[0, Length[gens]], "Verbose" -> False], "PIM"];
  ok = KeyExistsQ[$analyses, key] && $analyses[key]["detLogs"] === {} && $analyses[key]["rootLogs"] === {};
  ...
```

It runs `iPIM` on the residual as a probe and then reads the analysis back out of
`$analyses`. If that nested `iPIM` takes an early exit — `AnalyseSteps` `Throw`s
a certificate or a failure — then `$analyses[key]` is **never written**, and `ok`
is False. "Unverified" is then reported for a residual that may well be
residue-free: the probe simply did not complete.

That makes the answer depend on `$analyses` cache state, and `Part2` wraps each
run in `Block[{$analyses = <||>}, …]`, so it depends on what ran before.
Observed directly: on the same integrand and the same tower,
`VerifiedResidueFree[rem, T]` returns **True** when called eagerly and **False**
when evaluated lazily in the certificate's `&&` chain. With it forced True, every
other conjunct of the certificate is already satisfied on the outer run
(`OptionValue["Bounds"] === None`, `proved = True`, `typeE = False`,
`splittable = False`, `unitsComplete = True`).

### The plan

1. **Make `VerifiedResidueFree` answer the question it is asked.** Separate the
   three outcomes it currently collapses into one boolean: *residue-free*,
   *has residues*, *could not decide* (the probe threw, or the analysis is
   absent). A certificate must require the first; today it accepts only the first
   but reports the third as if it were the second. Concretely: have the probe
   return the analysis (or a status) rather than communicating through
   `$analyses`, so the result cannot depend on cache state or on a `Block` that
   reset it. `AnalyseSteps` already returns the Association — that is what the
   v0.239 refactor made possible, and this is its first consumer.
   *Reproduce with:* `ParallelMixed`Private`VerifiedResidueFree[{t, 0}, T]` for
   the tower `BuildTower[Exp[-x^2], x]`, called eagerly and then inside an `&&`.

2. **Re-validate and then delete the `q =!= None` guard.** Re-run the PMT stress
   corpus (`mixed/stress/run_stress.py --system mathilda`, 114 cases) and the
   review corpus with the guard lifted, and confirm T2, T10 and every other
   curve-free `T`/`M` case still answers rather than certifying. The guard should
   go only on that evidence, never on the two cases above alone — the whole point
   of the guard is that a certificate is unfalsifiable from inside.

3. **Then re-measure.** Expected: bucket A (39) and bucket B (11) close, since B
   is the partial mode declining to absorb exactly the sub-integrands of A. That
   is **247 → ~297**.

4. **Keep the soundness check honest while doing it.** Add the elementary
   curve-free integrands to a regression that asserts they are NEVER certified:
   `1/(x (Log[x]^2+1))`, `1/(x^4-1)`, `1/(x Log[x])`, `Exp[x]/(1+Exp[x])`,
   `Log[x]^2`, `x Exp[x]/(1+x)^2`. `tests/test_parallelmixedspecial.c` already
   has the shape for this in `test_certificate_withheld_honestly`; invert it once
   the certificate is expected to fire, so the test asserts *correct* certificates
   and *no* false ones rather than the current "withheld honestly".

---

## 3. Work item 2 — coverage (5 cases, bucket C)

Ids **87, 89, 90** (D trigonometric) and **172, 175** (P families) answer
`{"failed", "no solution with the special kernels"}` — the kernel search itself
does not close them. These are the only cases in the corpus where Mathilda's
*algorithm* falls short rather than its certificate, so they are the honest
measure of port fidelity: 5 of 312.

Start by diffing the verbose traces against Mathematica's for one of them
(`ParallelMixed`Private`IntegrateSurfaceSpecial[f, x, "Verbose" -> True]` prints
the Python's trace lines, so the three ports are line-comparable). The likely
suspects, in order: a missing `Ci`/`Si` conjugate-`Ei`-pair reduction in
`EiFromResidues` (S2), or the `Present` rewrite not recognising the pair.

---

## 4. Work item 3 — performance (2 cases, bucket D)

Ids **109, 110**:

```
1/((x^2 - 3x + 1) Sqrt[x^3 - x])        Python 45.0 s   Maxima 42.9 s   Mathilda >120 s
1/((x^2 - 2x - 1) Sqrt[x^3 - x])        Python 43.5 s   Maxima  ~3 s    Mathilda >120 s
```

Both need the **mod-p non-torsion certificate** of a third-kind residue divisor,
i.e. Jacobian arithmetic over GF(p). This is the gap already measured in
`MATHILDA_DIVERGENCES.md` §C: one group operation in `Jac(y^2 = f)` over GF(17)
costs Mathilda 10–20 ms against Python's 0.1 ms and Mathematica's 0.5 ms, and a
genus-4 certificate 512 s against Python's 0.8 s. The named fix is in the C core
and is unchanged by this stage:

> implement `Modulus -> p` for `PolynomialQuotient`, `PolynomialRemainder`,
> `PolynomialQuotientRemainder`, `PolynomialExtendedGCD` and `PolynomialMod` on
> FLINT's `nmod_poly` (already linked), so the straightforward Wolfram-language
> port runs as is

`PolynomialExtendedGCD[…, Modulus -> p]` and the `PolynomialQuotient` family were
fixed for correctness at v0.170 (A2/A3), but the `.m` still routes through
`ReduceModP`/`PolyQuoModP` — expression-level round trips — because the cost, not
the correctness, is what blocks it. Removing those workarounds in favour of the
now-correct builtins is the smallest useful step; it is worth measuring before
assuming a native path is needed.

Note this is worth **2 corpus cases** and is the least valuable item per unit of
work. It is listed third for that reason.

---

## 5. Work item 4 — branch fidelity (1 case, bucket E)

Id **34**, `asin(sqrt(x+1))/sqrt(x)`, returns `verified = False`. The residual is
**purely imaginary** (`0.0 - 1.29 I` at `x = 7/5`, `0.0 - 1.11 I` at `9/4`): the
real part is exact and the `ArcSin` branch in the reconstructed surface sits on
the other side of a cut from the integrand's.

**Part II alone fails the same integrand** (`ParallelIntegrateMixed` →
`{"failed", "verification failed"}`), so this is a pre-existing branch-fidelity
gap in the surface reconstruction, not something the stage introduced, and it
should be fixed there. The stage behaves correctly in reporting the answer
unverified rather than claiming it. Related prior art: the branch-resolved verify
gate of v0.175 (changelog §E.3), which pinned the sign numerically for A11/A34 of
Charlwood by accepting only exactly ±f at every real sample — the same treatment
extended to a complex sample point is the obvious first thing to try.

---

## 6. Expected arc

| after | PASS | what closed |
|---|---|---|
| today (v0.239) | 247 | — |
| item 1 (certificate) | ~297 | A 39 + B 11 |
| item 2 (coverage) | ~302 | C 5 |
| item 3 (GF(p) speed) | ~304 | D 2 |
| item 4 (branch) | **305** | E 1 |

305 is parity: the remaining 7 are the HONEST cases every port declines
(#59, #97, #98, #274–277 — strict mode withholding an uncertified special answer,
by design and identically in all four ports).

Item 1 is worth 50 cases and is the only one that needs a decision about
soundness rather than effort; items 2–4 are worth 5, 2 and 1.

---

## 7. How to re-measure

```bash
cd <research>/special
MATHILDA_BIN=/Users/user/Mathilda_dev/Mathilda/Mathilda \
  TIMEOUT=120 WORKERS=6 python3 stress_mathilda.py          # -> stress_mathilda.json
python3 stress_compare.py --ports=python,wl,maxima,mathilda # four-port table, exits 1 on any FAIL/WEAK
python3 stress_mathilda.py --only=C                         # one group
python3 stress_mathilda.py --ids=66,79,94                   # a selection
```

A copy of the runner and of the v0.239 results lives in `mixed/stress/`
(`stress_mathilda_special.py`, `stress_mathilda.json`,
`stress_report_mathilda.txt`) — note that directory is **untracked**, as the rest
of `mixed/` is.

Three measurement traps, each of which cost a debugging cycle and is now handled
in the runner but will bite any new harness:

- **`AbsoluteTime[]` has integer-second resolution** (divergence F7). The research
  runners time each case as an `AbsoluteTime[]` delta; here that reads `0.0` for
  anything under a second. Use `AbsoluteTiming`.
- **Two contexts** (divergence F8). Part II exports seven *public* symbols
  (`Tower`, `TowerD`, `ClassifyPrime`, `CertifyNonconstant`,
  `ParallelIntegrateMixed`, `BuildTower`, `Undecided`); everything else, and the
  whole special stage, is `ParallelMixed`Private``. Addressing either with the
  wrong context leaves the call **unevaluated**, which reads downstream as a
  wrong *answer*, not an error — it made all 41 group-A cases FAIL.
- **Do not compare summed times across ports.** 50 cases bail out early at the
  certificate step instead of doing the work, so Mathilda's total measures less
  work than the others'. Compare on a common subset. On the 247 cases Mathilda
  passes: Mathilda 41.2 s summed / 0.070 s median, Mathematica 56.0 / 0.110,
  Python 188.2 / 0.390, Maxima 569.7 / 1.930 — median per-case ratio 0.62×
  Mathematica, 0.16× Python, 0.04× Maxima. Speed is not the problem here;
  the certificate is.

## 8. Reference

- Paper: `<research>/special/rn-radicals-specfun.{tex,pdf}` (boxes S1–S10).
- Ports and their divergences: `<research>/special/PORTS.md`.
- Mathilda's own divergence log: `MATHILDA_DIVERGENCES.md`, section F for what
  this port exposed (F1–F6 fixed in v0.238, F7–F8 recorded).
- Session notes, including the two splice bugs the Part II merge hit and how they
  were caught: `tasks/parallelmixedspecial.md`.
- Changelog: `docs/spec/changelog/2026-09-28.md`, entries v0.238 and v0.239.
