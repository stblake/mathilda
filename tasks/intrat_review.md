# Review: the rational-function integrator (`src/calculus/intrat.c`, `intsimp.c`)

2026-11-27. Scope: `src/calculus/intrat.c` (~5000 lines), `src/calculus/intsimp.c`
(~2000), the `Method -> "BronsteinRational"` dispatch and the Automatic cascade in
`src/calculus/integrate.c`, checked against Bronstein, *Symbolic Integration I*,
ch. 2. Measured on the new benchmark `benchmarks/266-rational-integration` (389
cases, 15 families; see its `REPORT.md`).

**Line numbers** are those of the reviewed tree: HEAD 9338daf72 plus the main
tree's uncommitted patch (`../main_uncommitted.patch`). Where a finding was fixed
in this change set, the status says where the fix lives now.

## Outcome on benchmark 266 (Automatic)

First pass, before → after (389 cases):
- **Verified answers:** 367 → 389 of 389.
  - The 389th, `alg_root-011`, closes in 0.6 s with the merged Root-coefficient fix.
  - All 17 wrong answers (F1) and the 4 timeouts of F9 are gone.
- **Spurious discontinuities:** 4 → 0. All four were on the wrong answers.
- **Median size vs the reference:** nonsqfree 2.35× → 0.84×, multiplicity 1.79× →
  1.03×.

Second pass, before → after2 (423 cases; `czichowski` and `logtoreal` are new
families):
- **Verified:** 423 of 423, with 0 spurious discontinuities.
- **Median time vs the reference:** 0.01–0.34× in every family.
- **Median time per integral, before → after2:**
  - classical 4.0 → 1.7 ms;
  - nonsqfree 51 → 6.0 ms;
  - multiplicity 23 → 4.7 ms;
  - rootsum 65 → 4.7 ms;
  - sqfree 46 → 9.5 ms.
- **Worst case:** stress-019, an irreducible degree-20 denominator, went from 3.9 s
  in the first pass to 1.06 s (Czichowski); the reference takes 0.49 s.
- **Median size vs the reference:** ≤ 1.08× in every family except alg_root (1.85×,
  the Root-coefficient path).

## Second pass (2026-11-27): the four missing algorithms, and measured defaults

All four missing algorithms of Bronstein ch. 2 are now implemented natively in C, in
the new `src/calculus/intrat_alg.c`. They use FLINT `fmpq_poly` / `fmpz_poly` /
`fmpq_mat` / `qqbar` and the in-tree FLINT Gröbner engine, and the evaluator is
touched only to convert the input and the answer:
- **Quadratic Hermite reduction (§2.2):** `hermite_quadratic`. Partial fractions over
  the squarefree factorisation, then one ExtendedEuclidean per power.
- **Horowitz–Ostrogradsky (§2.3):** `hermite_horowitz`, one `fmpq_mat_solve` of size
  deg D.
- **Mack's linear Hermite (§2.2):** `hermite_mack`, native over Q.
- **Rothstein–Trager (§2.4):** `logpart_rt`. The resultant comes from deg D + 1
  integer evaluations and interpolation, skipping the one t at which A − tD′ drops
  degree. It is factored with `fmpz_poly_factor`, then gcd(D, A − aD′) is taken by
  Euclid over Q[t]/q for each irreducible factor q.
- **Czichowski (§2.6, Theorem 2.6.1):** `logpart_czichowski`. It computes the
  reduced lex basis of (D, A − tD′) with x > t on `gb_flint_buchberger`; then
  Q_i = content_x(P_i)/content_x(P_{i+1}) and S_i = pp_x(P_{i+1}).
- **General LogToReal (§2.8):** `intrat_alg_log_to_real`.
  - It finds the exact roots of every residue factor (`qqbar_roots_fmpz_poly`).
  - A real root a gives a Log[S(a, x)].
  - A root u + I v with v > 0 gives u Log[A² + B²] + v LogToAtan(A, B), where
    S(u + I v, x) = A + I B.
  - LogToAtan is Rioboo's continuous recursion, run exactly on qqbar polynomials
    (`ia_log_to_atan`).
  - This is the variant Bronstein recommends in practice (p. 68): compute the roots,
    then call LogToAtan. Theorem 2.8.4 guarantees gcd(A, B) = 1, so the recursion is
    well defined.
  - A coefficient beyond a quadratic surd is written as a `Root` object.

Each is selectable through `Method -> {"BronsteinRational", "Hermite" -> "Linear" |
"Quadratic" | "HorowitzOstrogradsky", "LogPart" -> "LazardRiobooTrager" |
"RothsteinTrager" | "Czichowski", "RealForm" -> "Rioboo" | "LogToReal" | None}`.
The same rules are accepted as trailing options of `Integrate`BronsteinRational`.
Plumbing:
- `intrat_method_suboption` / `_install` / `_restore` in `intrat.h`;
- `parse_method_option` in `integrate.c`;
- a setting is valid only for the top-level evaluation that installed it, so an
  abort cannot leak it into the next input.

Over Q every algorithm is native. Parameters, radicals and Root coefficients decline
to the Expr Mack reduction and LRT. That is where the number-field engine
(`src/poly/nfrat.c`, being built by the other agent) plugs in: `fmpq_poly` becomes
polynomials over Q(θ) in the same routines.

### Measurement: all 36 combinations on 423 cases (benchmark 266, `REPORT.md` "Algorithm grid")

The corpus gained two families for this measurement:
- `czichowski` (14 cases): high-degree squarefree denominators with few distinct
  residues;
- `logtoreal` (20 cases): quartics irreducible over Q that split over Q(√d) into
  complex quadratics, plus random cubics and quartics.

Each case runs all 36 configurations in one process, with 15 s per integral.
Verified means `PossibleZeroQ[D[F, x] − f]` and no spurious jump:

| configuration | verified | total time (406 common cases) | median |
|---|---|---|---|
| Linear · LRT · Rioboo (the previous default) | 423 | 17.43 s | 7.7 ms |
| Linear · RothsteinTrager · Rioboo | 423 | 14.30 s | 5.8 ms |
| Linear · Czichowski · Rioboo | 423 | 12.54 s | 5.9 ms |
| Quadratic · RT · Rioboo | 423 | 14.56 s | 5.4 ms |
| HorowitzOstrogradsky · RT · Rioboo | 423 | 14.02 s | 5.6 ms |
| any · any · LogToReal (Root objects) | 409–422 | 15.7–24.2 s | 7.4–11.6 ms |
| any · any · None | 397–398 (22–23 JUMP: complex logs, expected) | 10.9–19.4 s | 4.8–8.5 ms |

Per family (from `REPORT.md`):
- **Hermite.** The three variants agree on every answer and are within noise of each
  other in every family.
- **Log part.**
  - RothsteinTrager is the fastest in 13 of 17 families. Examples: largecoef 3.1 s
    against LRT 4.4 s and Czichowski 5.4 s; nonsqfree 0.35 s against 0.52 s.
  - Czichowski is 4–5× faster than both alternatives where the residue factors are
    irreducible of high degree: rootsum 0.44 s against RT 1.99 s and LRT 1.76 s;
    stress 1.08 s against 5.05 s and 5.19 s.
- **Real form.**
  - Rioboo (the shape dispatch) never lost a case.
  - The general LogToReal with Root objects verifies, but its answers are 1.4–4×
    larger. It is also expensive on high-degree factors: sqfree takes 82 s.
  - Automatic's experimental "general LogToReal in radicals when smaller" never
    produced a smaller answer on the corpus.

**Chosen defaults (Automatic):**
- `"Hermite" -> "Linear"` (Mack), natively over Q.
- `"LogPart"`:
  - RothsteinTrager while every irreducible residue factor has degree ≤ 8;
  - else Czichowski;
  - else LRT, which also covers parameters, radicals and Root coefficients.
- `"RealForm" -> "Rioboo"`.
- The alternatives are not run as a cross-check by default. They agree on all 423
  cases, and the cheapest independent check, `PossibleZeroQ`, is what the benchmark
  already runs.

### Findings of the second pass

- **S1. The general LogToReal needs size caps before any exact-root arithmetic.**
  qqbar operations cannot be interrupted by `TimeConstrained`. With 20- to 60-digit
  coefficients, Re and Im of a cubic's roots produced `Root` objects with
  thousand-digit minimal polynomials. `PossibleZeroQ` could not decide those, and
  one runaway killed the whole grid process.
  - Caps: residue factor degree ≤ 6 and coefficients ≤ 256 bits with Root objects;
    degree ≤ 4 in radicals.
  - Rothstein–Trager declines when a factor has degree > 8: the degree-20 stress case
    took more than 15 s in the gcd over Q(a).
- **S2. Process-order warm-up biases in-process A/B timing.** The first configuration
  in a process pays the lazy initialisation, which made the "Linear" Hermite look 2×
  slower even on squarefree inputs. The grid tables compare configurations that
  never run first.
- **S3. Pre-existing leak (not in this code).** `PossibleZeroQ[D[F, x] − f]` on
  `F = RootSum[1 + 3 #1 − 31 #1^3 &, Log[4 + 9 x + 31 #1 − 62 #1^2] #1 &]` leaks 10
  blocks (480 bytes) per call in the previous binary too. The leak is in the
  RootSum derivative / zero-test path. The new code paths are leak-free
  (`leaks --atExit`, every algorithm).
- **S4. `RealForm -> None` now prints in the reference's slot form**
  (`RootSum[... &, Log[...] #1 &]`), not `Function[tt$, ...]`.

## What is implemented (against Bronstein ch. 2)

| Bronstein ch. 2 | in Mathilda (after the 2026-11-27 second pass) | notes |
|---|---|---|
| §2.2 Hermite reduction, quadratic | **yes**, native: `hermite_quadratic` (`intrat_alg.c`), `"Hermite" -> "Quadratic"` | partial fractions over the squarefree factorisation, then one reduction per power |
| §2.2 Hermite reduction, linear (Mack) | **yes**, native over Q: `hermite_mack` (`intrat_alg.c`); the Expr port `intrat_hermite_reduce` serves other coefficient domains | the Automatic default |
| §2.3 Horowitz–Ostrogradsky | **yes**, native: `hermite_horowitz` (`fmpq_mat_solve`), `"Hermite" -> "HorowitzOstrogradsky"` | `hermiteRat` in `goursat/RatInt.m` still duplicates it in `.m` |
| §2.4 Rothstein–Trager | **yes**, native: `logpart_rt` (resultant by evaluation/interpolation, `fmpz_poly_factor`, Euclid over Q[t]/q), `"LogPart" -> "RothsteinTrager"` | the Automatic default over Q while every irreducible factor has degree <= 8 |
| §2.5 Lazard–Rioboo–Trager | **yes**, `intrat_log_part_core` on the native subresultant chain (`intrat_subres.c`) | parameters, radicals, Root coefficients, factors of degree > 8 |
| §2.6 Czichowski (Gröbner) | **yes**, native: `logpart_czichowski` on `gb_flint_buchberger` (lex, x > t), `"LogPart" -> "Czichowski"` | fastest on long chains with few residues (stress-018 0.17 s vs LRT 1.1 s, RT 10 s), not the default (see the grid) |
| §2.7 Newton–Leibniz–Bernoulli (full partial fractions) | implicitly: `intrat_integrate_summands` `Apart`s h over Q before the log part | (F10) |
| §2.8 LogToAtan (Rioboo, continuous) | **yes**: `intrat_log_to_atan` (Expr) and `ia_log_to_atan` (exact, qqbar polynomials) | |
| §2.8 LogToReal | **yes**, general: `intrat_alg_log_to_real` -- exact roots of every residue factor (qqbar), real root a: a Log[S(a, x)], u + I v: u Log[A^2 + B^2] + v LogToAtan(A, B); `"RealForm" -> "LogToReal"` | the shape dispatch `logtoreal_dispatch` remains as "Rioboo"; Automatic adds the general one in radicals when it is smaller |
| Verification of the answer | none | (F2) |

## Findings, ranked by impact

Status: **FIXED** (in this change set), **OPEN**, **OWNED** (another agent's area:
the Root-coefficient log part and `intrat_reduce_root_coefficients`).

### Correctness

**F1. A BigInt discriminant was read as sign-unknown, giving a wrong answer.** — FIXED
`intsimp_sign_pos_assumption` (`intsimp.c:115-198`) handled `Integer` and a
`Rational` of two machine integers only; `EXPR_BIGINT` and a `Rational` with a BigInt
part fell through to 0 ("undecided"). `logtoreal_quadratic_core` (`intrat.c:3255`)
then took the complex-pair branch for a *positive* discriminant, with an imaginary
`Sqrt[-disc]`, and returned a single wrong `Log` of a quadratic with huge nested
radicals. This happens whenever the resultant's coefficients pass 2^63, which
Hermite reduction produces readily. It caused all 17 FAILs of the baseline: 5 of 30
in nonsqfree and 12 of 25 in largecoef. Repro:
`(9814156935056 - 2671246202843 x)/(1252991435641 (x^2 + 4 x - 1))`.
The same function also could not sign a parameter-free algebraic constant such as
`3 - 2 Sqrt[2]`: the Plus rule returns 0 for mixed-sign terms. Such inputs got a
formally correct but complex-valued `Sqrt[-3 - 2 Sqrt[2]] ArcTan[...]`.
*Fix:* Integer/BigInt/Rational signs via `mpz_sgn`, then a certified ball sign
(`flint_cert_parts_sign`) for parameter-free exact constants. A non-real constant
answers 0, so it keeps the formal branch-agnostic form.

**F2. Log terms silently dropped when an S_i cannot be found.** — FIXED
In `intrat_log_part_core` (`intrat.c:2118`, `:2128`, `:2155`), three failures only
`continue`: no PRS element of x-degree i, `SquareFree` of lc(S) failing, and
`intrat_primitive_part_mod` returning NULL. The output loops (`:2169`, `:2219`) then
skip that Q_i, so its entire RootSum vanishes from the antiderivative. Bronstein's
IntRationalLogPart guarantees an S_i for every nontrivial Q_i. *Fix:* the core now
declines (NULL) when any Q_i of positive degree has no S_i. Nothing in `intrat.c`
verifies D[F] = f, so this was a silent wrong answer.

**F3. `intrat_exquo` returned an inexact quotient.** — FIXED
At `intrat.c:166-176` a nonzero remainder printed a raw `fprintf(stderr,
"Integrate`exquo::remainder ...")`, which bypasses `Quiet`/`Check`. It then returned
the quotient anyway. Hermite (`:569`, `:587`, `:597`, `:620`), the S_i correction
(`:2145`) and the partial log part (`:2075`) carried on with wrong data. *Fix:*
return NULL, which every caller already treats as "cannot integrate".

**F4. LogToReal over a non-real coefficient field.** — FIXED (as part of F1)
The guard in `try_lrt_close` (`:4312`) looks only for literal `Complex[]` leaves, so
`(-1)^(1/3)` or a complex `Root` passed it, and `split_re_im` treated it as real. With
F1's certified sign, a non-real discriminant no longer reads as "positive". Checked:
`1/(x^2 + (-1)^(1/3))`, `1/(x^2 + x + (-1)^(1/3))` and `1/(x^3 + (-1)^(1/3))` all
verify. The literal-`Complex` gate itself is still narrow (OPEN, low).

**F5. `LogToAtan` never checks that its ArcTan argument is a polynomial.** — OPEN (medium)
`v = canonic((AD+BC)/g)` at `:2569-2575` is assumed to be a polynomial (Bronstein
§2.8, p. 63). When `PolynomialExtendedGCD` misses a common factor that holds only
through a symbolic `Sqrt` or a parameter, v is rational. ArcTan[v] then jumps by π
at v's real poles: a spurious discontinuity, with the derivative still correct. No
case in the corpus triggers it. Fix: decline, or keep the RootSum, when v is not
polynomial in x.

**F6. Log merging and content stripping on complex kernels.** — OPEN (medium)
`intsimp_log_to_arctanh` (`intsimp.c:594`) merges c Log[A] + c Log[B] into
c Log[A B], and `intsimp_strip_log_constants` drops `Log[c]`. Both are exact up to
piecewise constants for real A, B. On the I-coefficient path (Log[x − α] with α
complex) they add 2πi jumps at regular real points. Gate both on real coefficients.

**F7. Parameter assumptions are implicit.** — OPEN (low-medium)
`intsimp_sign_pos_assumption` treats every symbol as positive. Leading-coefficient
divisions (`intrat_polygcd_monic` `:297`, `intrat_extract_constants` `:1088`) assume
generic parameter values. The forms produced (two-root, ArcTan, biquadratic,
nth-root) are formal identities for any branch of the square root. So the
derivative is right for every parameter value; only real-valuedness is assumed. No
condition or message tells the user this. The special-value family (`param_special`)
passes 20/20.

**F8. ArcTanh oddness and ArcCoth probes.** — OPEN (low)
`intsimp_normalize_inverse_trig_signs` (`intsimp.c:835`) uses ArcTanh[−z] =
−ArcTanh[z], which fails by iπ on the cut. `it_arg_real_outside_unit`
(`intsimp.c:1270`) probes only x > 0.

### Speed

**F9. The "subresultant" chain was an exponential pseudo-remainder sequence.** — FIXED (largest speed item)
`intrat_log_part_core` took its chain from the `SubresultantPolynomialRemainders`
builtin (`poly.c:3171-3233`). That builtin is a plain pseudo-remainder sequence with
no content removal and no subresultant scaling, run on Expr trees, so its
t-coefficients grow exponentially (Bronstein §1.5). For Czichowski's degree-14
example (corpus `bronstein-009`) the chain had a leaf count of 71,924. It took 19.2 s
of the 18.9 s integral, while `Resultant` took 2 ms. The baseline's 4 stress
timeouts (`1/(x^30+1)` and the irreducible degree-12/16/20 denominators) are the
same cost.
*Fix:* new `src/calculus/intrat_subres.c` computes the genuine subresultant PRS
(Bronstein's SubResultant, Brown–Collins γ/β scaling, every division exact) on dense
FLINT arrays: x is the main variable, each coefficient an `fmpz_poly` in t. It
returns the same chain elements up to constants in Q(t), so the downstream reads are
unchanged. Scope: exact rational coefficients; parameters, radicals, Root objects and
the transcendental callers keep the generic path. Results: Czichowski 18.9 s →
0.025 s; `1/(x^30+1)` >30 s → 0.11 s; the irreducible degree-12 case >30 s → 0.6 s.
*Remaining:* the generic path for parameter or algebraic coefficients still uses the
exponential builtin. The builtin should become a real subresultant chain: Ducos, or
the same algorithm over `fmpz_mpoly` / `gr_poly`. It lives in `poly/`, owned by the
Together/Factor agent this week. `find_prs_at_degree` (`:1254`) re-expands every
chain element for each i, which is quadratic but harmless now.

**F10. Full Apart before LRT, then ~5 factorisations per piece.** — OPEN (medium-high)
`intrat_integrate_summands` (`:4330`) `Apart`s h over Q, which factors the
denominator: the Newton–Leibniz–Bernoulli approach of §2.7 that LRT (§2.5) exists to
avoid. Each piece then gets:
- two `FactorSquareFree` calls in `ExtractConstants` (`:1076-1079`);
- one in `derivative_recognition` (`:786`), already done on the whole denominator;
- the resultant's squarefree split (`:2025`), and one per lc(S_i) (`:2126`);
- `Factor` of Q_i in LogToReal (`:4049`), and again in the linear closer (`:4157`).

The per-factor split does give the real forms the reference prints, so the fix is to
share one factorisation rather than drop the split.

**F11. Full `Simplify` inside the S_i normalisation.** — OPEN (medium)
`intrat_primitive_part_mod` (`:518`) calls `Simplify` once per Q_i whenever a radical
is present. The normalisation is optional (Bronstein p. 51). Reduce through
`flint_algebraic_field_normalize` instead.

**F12. Redundant lc(S_i) correction.** — OPEN (low)
After `s = intrat_primitive(s_raw, x)` (`:2119`) has removed the whole t-content,
the gcd(A_j, Q_i) correction (`:2122-2148`) cannot fire, yet it still costs a
`SquareFree` and q gcds. Keep one of the two.

**F13. Hermite reduction runs entirely through the evaluator.** — FIXED over Q (second pass: native `hermite_mack` / `hermite_quadratic` / `hermite_horowitz`); OPEN for parameters and algebraic coefficients
Per level it does two gcds, one xgcd, three exact quotients and two Together/Cancel,
and D[dbar] is computed twice (`:577`, `:595`). It recomputes the squarefree data
that `FactorSquareFree` (`:786`) and `intrat_squarefree_list` already have. Over Q,
`fmpq_poly_gcd/xgcd/divrem/derivative` from one `fmpq_poly_factor_squarefree` would
make it native. On the corpus Hermite is not the bottleneck: multiplicity cases are
25 ms median, 6.7× faster than the reference.

**F14. `intsimp_log_to_arctanh` is O(n²) with evaluator zero tests per pair.** — OPEN (low-medium)
At `intsimp.c:590`. Bucket the terms by coefficient first.

**F15. Four Together passes on entry.** — OPEN (low)
`integrate.c:247`, `intrat.c:4412`, `:2254`, `:1712`.

### Output compactness

**F16. The Hermite rational part was split over its denominator.** — FIXED
The driver passed g, a sum of per-level `b/dbar` fractions, through
`intsimp_distribute_plus`, which is `expr_expand`. That split each numerator over
its expanded denominator: `63936/D - 165024 x/D - 2592 x^2/D`. *Fix:* the new
`intsimp_compact_rational` (`intsimp.c`, end) picks the smallest LeafCount of three
forms: one fraction over the squarefree-factored denominator, Apart over Q, or g as
given. g is now kept out of the log-part passes. Examples: `1/(x^2+1)^8` 195 → 52
leaves; `1/(x^3+x+1)^4` 313 → 92; Bronstein's Hermite example 114 → 64 (the
reference prints 68).

**F17. Log arguments kept their rational content.** — FIXED
`intsimp_strip_log_constants` only rewrote `Log[Times[c, p]]`, so
`Log[4 - 4x + 4x^2 - 4x^3 + 4x^4]`, `Log[25 + 25 x + 25 x^2]` and
`Log[-148204931958652908 + 49401643986217636 x + ...]` kept their content. *Fix:* a
`Log[Plus[...]]` is divided by the rational content of its coefficients (BigInt
safe).

**F18. When LogToReal fails, the LRT result is thrown away.** — OPEN (medium)
For a Q_i of degree ≥ 3 that is not one of the special shapes, `try_lrt_close`
(`:4321`) discards all pairs and returns `intrat_naive_log_part`: a RootSum over the
whole piece denominator, with residues a(α)/d′(α) and Log[x − α]. The minimal
Σ RootSum[Q_i, t Log[S_i]] is already built by `intrat_root_sum_pairs` (`:2370`) on
the transcendental path. The fallback is also all-or-nothing: one bad Q_i discards
the real forms of the other pairs. Fix: fall back per pair. On the corpus this
mostly costs nothing, because the piece denominators are irreducible so
deg Q_i = deg d, but it is the general case.

**F19. LogToReal lacks the general algorithm.** — FIXED (second pass: `intrat_alg_log_to_real`, see the table above and the grid section of benchmark 266)
Bronstein's general LogToReal (§2.8, p. 69) substitutes u + iv and solves
P(u, v) = Q(u, v) = 0. It is missing. The palindromic r = 1 case duplicates
`expand_palindromic_quartic_real` (`:1351`). The comment at `:3836` claims it catches
`1/(x^5+1)`, but the r = 1 restriction (`:3872`) makes it bail. The sparse nth-root
case (`:3417`) is justified: it gives the reference's `Cos[kπ/n]` form, e.g.
`1/(x^16+1)`, 577 leaves against the reference's 417.

**F20. Equal-coefficient logs merged blindly into an expanded product.** — OPEN (low)
`intsimp.c:594`; gate it on the merged form being smaller.

### Dead code, duplication, stale text

- **F21.** `intrat_nthroot_q` (`:1143`, `MATHILDA_MAYBE_UNUSED`) is dead.
  `intrat_has_complex` (`:2267`) and `intrat_contains_complex` (`:3214`) are
  identical. `hermiteRat`/`algNorm` in `goursat/RatInt.m` duplicate Hermite and the
  norm log part.
- **F22.** Leak: `:1668` and `:1670`, `nlp_sub_body(..., eval_and_free(nlp_neg(...)))`
  copies its value argument, so the negated root leaks.
- **F23.** The options `PFD`, `LogToArcTan`, `Radicals` and `Extension` are accepted
  and silently ignored (`:4564-4610`).
- **F24.** Stale comments:
  - ~63 references to line numbers of `IntegrateRational.m`, which is not in the
    repo;
  - the `intrat_canonic` comment about Extension;
  - `:4209` says "degree >= 5" where it means ≥ 3;
  - `:2654` says "fall through to RootSum" where it means the naive part;
  - `:4182-4198` claims a whole-h fallback that does not exist (the code returns
    NULL).

### Dispatch and cascade (`integrate.c`)

**F25.** The rational stage is reached early in the cascade: after Undefined,
UnitStep and the elliptic-trig check (`:4839-4849`), and after
`integrate_split_constant`. `Method -> "BronsteinRational"` calls `try_rational`
directly (`:5062`), then `intsimp_finalize` (`:5171`). On the whole corpus Automatic
and BronsteinRational return identical answers, so nothing earlier intercepts a
rational integrand. When the rational stage declines, the cascade continues through
LinearRadicals, Linearity (which re-runs the rational stage per term),
DerivativeDivides and Risch; only ParallelMixedTower skips rational integrands
(`:1767`). With F2/F3 now declining instead of answering wrongly, a decline costs the
whole cascade. Making the rational stage final for rational integrands, with the
naive RootSum as the terminal fallback, would bound that. OPEN (low-medium).

### Owned by the Root-coefficient agent (brief)

- **F26.** `intrat_reduce_root_coefficients` runs three times on every output
  (`:4293`, `:4512`, `:4537`).
- **F27.** The repro `r = Root[837/4 + #1^3 &, 1]`, ρ in x (corpus `alg_root-011`)
  times out at 30 s under both Automatic and BronsteinRational.

## A general mechanism for algebraic coefficients

The merged `intrat_log_part_core_root` (one Root coefficient, number-field
arithmetic through `PolynomialExtendedGCD`/`QuotientRemainder` with `Extension`)
fixes the Root hang. It is still one coefficient kind with its own path, next to
`intrat_reduce_root_coefficients`, the `Sqrt`/`I` handling and the Simplify in F11.
The general form is the log part over one number field K = Q(θ):
- **The field.** Take a primitive element of all the integrand's algebraic
  constants, whether `Sqrt`, `I`, nested radicals or `Root` (`flint_qqbar` already
  computes minimal polynomials).
- **The subresultant chain** of `intrat_subres.c`, with `fmpz_poly` (Z[t]) replaced
  by K[t]: `gr_poly` over `gr_ctx_init_nf`, the layer `flint_numberfield_gcd`
  already uses. Its divisions are exact in K[t] just as they are in Z[t].
- **S_i normalisation and LogToAtan** in the same context.

This would subsume `intrat_log_part_core_root`, the Root reduction passes (F26) and
the Simplify in F11, and would remove the last exponential path (F9's generic chain)
for every algebraic input. Parameters would need `fmpz_mpoly_q` coefficients in
place of K, which is a separate step. Recommended as the next piece of work; not
done here.


**Update (nfrat merged).** The number-field engine `src/poly/nfrat.c` is now the
substrate for this. It covers K = Q(I, Sqrt[n], n^(1/k), (-1)^(p/q), Root), any
number of generators, with a certified gcd over K. Today it exports only
Together / Cancel, and rat.c routes the Expr path of `intrat.c` (Mack's Hermite
and the LRT S_i normalisation, via `intrat_canonic`) through it. That made the
algebraic families faster: alg_root median 30 ms → 18 ms, alg_qi 14 → 8 ms. The
native routines of `intrat_alg.c` (all three Hermite variants, RT, Czichowski,
LogToReal) need nfrat to export its K[x] layer: conversion from Expr, divrem, gcd
and xgcd on its fmpq_mpoly-mod-theta polynomials. With that layer:
- `fmpq_poly` becomes K[x];
- Q[t]/q becomes K[t]/q;
- the qqbar roots of LogToReal come from the norm of each residue factor.

No second mechanism should be added; the export is the next step, and it belongs to
nfrat's owner.

## Where native FLINT arithmetic would pay most (after this change set)

1. **The generic subresultant chain** for parameter and algebraic coefficients
   (`poly.c:3171`), via `fmpq_mpoly` or `gr_poly`. It is the last exponential path.
2. **Hermite reduction over Q** (`:545-660`) with `fmpq_poly`, sharing one
   squarefree factorisation with the log part.
3. **S_i normalisation modulo Q_i** (`:2101-2160`, `:507-538`) with
   `fmpq_poly_xgcd`/`rem`, or `nf_elem` when Q_i is irreducible; this also removes
   the `Simplify`.
4. **The `LogToAtan` recursion** (`:2486-2605`) with `fmpq_poly_xgcd` over Q, or
   `gr_poly` over `gr_ctx_init_nf`.
5. **One shared factorisation** in place of the ~5 per Apart piece (F10).
