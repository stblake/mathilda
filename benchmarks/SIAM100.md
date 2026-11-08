# The SIAM 100-Digit Challenge (experiments 166–195)

Source: F. Bornemann, D. Laurie, S. Wagon, J. Waldvogel, *The SIAM 100-Digit
Challenge: A Study in High-Accuracy Numerical Computing* (SIAM, 2004). The book
covers the ten problems Trefethen set in *SIAM News* (Jan/Feb 2002), and its
Appendix D adds twenty "More Problems".

Each problem's answer is **one real number**. Mathilda is measured against four
other systems:
- **wolframscript**, which runs the same `.m` file;
- **Python**: scipy at machine precision, mpmath for the 50-digit row;
- **Maxima**;
- **FriCAS**.

Every problem has two tiers:
- **machine**: IEEE doubles, labelled `<id>.machine`;
- **wp50**: a 50-digit answer from arbitrary-precision arithmetic, labelled `<id>.wp50`.

There is also a reduced-instance row (`<id>.machine n=…`) where a double cannot
carry the full problem.

**Target.** On every case, Mathilda must:
- reach at least as many correct digits as the best other system; and
- run in under 1.5× the time of the fastest other *correct* system.

`siam100_table.py` scores this; its header gives the exact rules.

## Methodology

- **Idiomatic statements only.** Each system is asked the problem in its most
  direct built-in form: `NIntegrate`, `NMinimize`, `NDSolve`, `FindRoot`,
  `NSum`, `LinearSolve`, `Eigenvalues`… in WL; the scipy/mpmath equivalents in
  Python; QUADPACK, `find_root`, `rk`, `lbfgs`, `bfloat` in Maxima;
  `romberg`/`simpson`, `solve`/`complexSolve` and the `Float` domain in FriCAS.
  - Never the book's specialised algorithm: no contour trapezoid sums, AGM
    closed forms or hand-written interval Newton, unless that IS the direct
    statement.
  - The question is how good each system's *library* is.
  - The same mathematical formulation is used in every system. Only the
    spelling changes.
- **References.** Every true value is known to ≥ 55 digits. It is established
  either by two independent methods that agree, or (for Appendix D problems
  with no published value) by `consensus`: two different systems or algorithms
  at higher precision agree on the digits recorded, and no more digits than
  that.
  - The book prints only the first and last 25 of its 10,000 digits, so the
    50-digit references are computed here.
  - Provenance is on each line of the experiment's `truth.tsv`.
- **Correct.** At least `need` correct digits:
  - typically 10 at machine precision, matching the challenge's own scoring;
  - typically 48 for wp50;
  - lower only where the problem itself limits a double, for example C2,
    where 10 digits need roughly 40 digits of input.
- **Overfitting guard.** There are 10 generalised problems per problem in
  holdout A (`SIAM100_HOLDOUT_A.md`, experiments 196–225), and 10 more in a
  sealed holdout B (`SIAM100_HOLDOUT_B.md`), which is run once on the final code.

## File layout per experiment `NNN-siam100-<id>-<slug>/`

| file | system | helper loaded first |
|---|---|---|
| `<slug>.m` | Mathilda **and** wolframscript (same file) | `Get["../harness.m"]; Get["../siam100_common.m"];` |
| `<slug>.py` | Python 3.11: scipy (machine) and mpmath (wp50) | `from siam100_common import vrun, skip, require` |
| `<slug>.mac` | Maxima | `load("../harness.mac")$ load("../siam100_common.mac")$` |
| `<slug>.input` | FriCAS | `)read ../siam100_common.input )quiet` |
| `truth.tsv` | the true answers | label TAB value TAB need TAB provenance |

- **Case wrappers:** `vrun[label, expr, d(, reps(, cap))]` (WL), `vrun(label,
  fn, d, reps=1, cap=120, dps=None)` (Py), `vrun(label, lambda([], …), d)`
  (Maxima), and `t0 := nowms(); v := …; vcaseD/vcaseF(label, v, nowms() - t0,
  d)` (FriCAS).
- **Labels:** every half emits the **same labels**.
- **Unsupported cases:** a case that a system cannot state with its built-ins
  is `skip(label, "Head")`. It must never be a hand-written algorithm standing
  in for a missing built-in.
- **Output:** cases print `BENCH`, `CHECK` and `VALUE` lines. `VALUE` carries
  every digit and is what gets scored.
- **Runner:** `run_all.py` runs the `.m` file in both Mathilda and
  wolframscript, the `.py` file under `HPC_PYTHON`, the `.mac` file in Maxima,
  and the `.input` file in FriCAS (on stdin).

## The thirty problems

### Challenge problems (Trefethen, *SIAM News* 2002)

| exp | id | problem | answer begins |
|---|---|---|---|
| 166 | c1 | lim_{ε→0} ∫_ε^1 x⁻¹ cos(x⁻¹ log x) dx | 0.32336 74316 |
| 167 | c2 | A photon moving at speed 1 in the x-y plane starts at t=0 at (1/2, 1/10) heading due east. Around every integer lattice point (i,j) there is a circular mirror of radius 1/3. How far from (0,0) is the photon at t=10? | 0.99526 29194 |
| 168 | c3 | The infinite matrix A with a₁₁=1, a₁₂=1/2, a₂₁=1/3, a₁₃=1/4, a₂₂=1/5, a₃₁=1/6, … (a_ij = 1/((i+j−1)(i+j−2)/2 + j)) is a bounded operator on ℓ². What is ‖A‖? | 1.27422 41528 |
| 169 | c4 | Global minimum of e^{sin(50x)} + sin(60e^y) + sin(70 sin x) + sin(sin(80y)) − sin(10(x+y)) + (x²+y²)/4 | −3.30686 86474 |
| 170 | c5 | f(z) = 1/Γ(z); p(z) is the cubic that best approximates f on the unit disk in the sup norm. What is ‖f − p‖∞? | 0.21433 52345 |
| 171 | c6 | A flea starts at (0,0) on the 2-D integer lattice and does a biased random walk: N or S with probability 1/4 each, E with 1/4+ε, W with 1/4−ε. The probability that it ever returns to (0,0) is 1/2. What is ε? | 0.06191 39544 |
| 172 | c7 | A is the 20000×20000 matrix with the primes 2, 3, 5, …, 224737 on the diagonal and 1 at every a_ij with \|i−j\| = 1, 2, 4, 8, …, 16384. What is (A⁻¹)₁₁? | 0.72507 83462 |
| 173 | c8 | The square plate [−1,1]² is at u=0. At t=0 one side is raised to u=5 and held there, the other three stay at 0, and u_t = Δu. When does the centre reach u=1? | 0.42401 13870 |
| 174 | c9 | I(α) = ∫₀² (2 + sin(10α)) x^α sin(α/(2−x)) dx. For which α ∈ [0,5] is I(α) maximal? | 0.78593 36743 |
| 175 | c10 | A particle at the centre of a 10×1 rectangle does Brownian motion until it hits the boundary. What is the probability that it hits one of the ends rather than one of the sides? | 3.83758 79792 e−7 |

### Appendix D, "More Problems"

| exp | id | problem |
|---|---|---|
| 176 | d1 | Σ 1/n over the positive integers n whose decimal representation does not contain the substring "42". |
| 177 | d2 | Σ_{n≥1} 1/f(n), where f(1)=1, f(2)=2, and for n>2, f(n) = n·f(d) with d the number of base-2 digits of n. |
| 178 | d3 | E_N is the minimal value of Σ_{j<k} \|x_j − x_k\|⁻¹ over N points on the unit sphere. What is E₁₀₀? |
| 179 | d4 | Riemann's R(x) = Σ_{k≥1} μ(k)/k · li(x^{1/k}) (li as a principal value). What is the largest positive zero of R? |
| 180 | d5 | A is the 48×48 Toeplitz matrix with −1 on the first subdiagonal, +1 on the main diagonal and the first three superdiagonals, 0 elsewhere. What is min_p ‖p(A)‖₂ over monic polynomials p of degree 8? |
| 181 | d6 | ∫_{−1}^{1} exp(x + sin e^{e^{x+1/3}}) dx, i.e. exp(x + sin(exp(exp(x + 1/3)))) (the book's typeset tower read as two exponentials) |
| 182 | d7 | f(x,y) = e^{−(y+x³)²}, g(x,y) = y²/32 + e^{sin y}. What is the area of the region of the plane where f > g? |
| 183 | d8 | c_m(ℝᴺ)² = (2π)^{−N} ∫_{ℝᴺ} (Σ_{\|k\|≤m} x^{2k})⁻¹ dx, with k a multi-index. For example c₃(ℝ¹) = 0.5. What is c₁₀(ℝ¹⁰)? |
| 184 | d9 | A particle starts at the top vertex of a triangular lattice with 30 points on each side and takes 60 random steps (uniform over its lattice neighbours). What is the probability that it ends in the bottom row? |
| 185 | d10 | x₀ = x₁ = 1, x_{n+1} = 2x_n ± x_{n−1}, with each sign chosen independently with probability 1/2. To what does \|x_n\|^{1/n} converge almost surely? |
| 186 | d11 | Six unit masses are fixed at (2,−1), (2,0), (2,1), (3,−1), (3,0), (3,1). A unit mass starts at (0,0) with speed 1 in direction θ and moves under the inverse-square attraction r⁻² of the six. What is the shortest time in which it can reach (4,1)? |
| 187 | d12 | T = ½(ẋ² + ẏ²), U = y + (ε⁻²/2)(1 + αx²)(x² + y² − 1)². The particle starts at (0,1) with velocity (1,1). For which α does it first hit y=0 at time 10, in the limit ε→0? |
| 188 | d13 | u=(x,y,z) starts at (0,0,z₀), z₀≥0, with ẋ = −x + 10y + ‖u‖(−0.7y − 0.03z), ẏ = −y + 10z + ‖u‖(0.7y − 0.1z), ż = −z + ‖u‖(0.03x + 0.1y). If ‖u(50)‖ = 1, what is z₀? |
| 189 | d14 | −Δu = exp(α‖x‖²) on the regular pentagon inscribed in the unit circle, with u=0 on four sides and a Neumann condition on the fifth. For which α does the integral of u along the Neumann side equal e^α? |
| 190 | d15 | At what time t_∞ does the solution of u_t = Δu + e^u on the 3×3 square, with zero boundary and zero initial data, blow up? |
| 191 | d16 | φ₂ is the Daubechies scaling function (support [0,3]). The heat equation u_t = u_xx on [0,3] is solved with u(x,0)=φ₂(x) and u(0)=u(3)=0. At what time does the length of the curve u(·,t) in the x-y plane fall to 5.4? |
| 192 | d17 | u is an eigenfunction for the third Dirichlet eigenvalue of the Laplacian on the L-shaped domain made of three unit squares. What is the length of its zero-level set? |
| 193 | d18 | The smallest eigenvalue of −Δ with zero boundary conditions on the Koch snowflake built on an equilateral triangle of side 1. |
| 194 | d19 | −div(c(x) grad u) = 1 on the unit square with zero Dirichlet data. What is the infimum of ∫u over c taking the value 1 on half the area and 100 on the rest? |
| 195 | d20 | h is the solution of exp(z) = h(h(z)) that is analytic near 0 and increasing on the reals. What is h(0)? (The radius-of-convergence question is open.) |

**Infeasible tiers.** Where no system can attempt a tier with its built-ins
(for example a 50-digit PDE eigenvalue on a fractal domain), the case is
written anyway and emits `SKIP` or a timeout in every system. It is listed
here as **N/A**. A tier is never dropped just because Mathilda fails it.

## Status

### Baseline, v0.6377 (2026-11-07), before any fix

Full table: `siam100_results/baseline_v0.6377.md`. Raw run: siam100_table.py
`--run --only 166-195`, serial, but concurrent with holdout-A authoring agents. So:
- timings are approximate;
- several wolframscript cells are "—" because of licence contention (c5, c6, c7,
  c10 machine), and will be re-run on a quiet machine.

Over 76 cases, Mathilda is on target in 48. Correct cases per system:

| system | correct | inexact | wrong | none |
|---|---|---|---|---|
| Mathilda | 25 | 9 | 8 | 34 |
| wolframscript | 21 | 10 | 5 | 40 |
| Python (scipy/mpmath) | 33 | 6 | 7 | 30 |
| Maxima | 20 | 2 | 1 | 53 |
| FriCAS | 14 | 0 | 2 | 60 |

Off target, grouped by cause (the work queue; see tasks/todo.md):
- **Wrong answer:**
  - c1.wp50: `NIntegrate` at WP 60 gives −1.19;
  - d19: constrained `FindMinimum` makes no progress.
- **Silently ignored option:**
  - d12, d13: `WhenEvent` in `NDSolve`;
  - d15: `WorkingPrecision` on a PDE.
- **MPFR tier absent or timeout:**
  - c4: `NMinimize` returns only machine digits;
  - c7: sparse `LinearSolve`;
  - c9: `FindMaximum` of `NIntegrate`;
  - c10: `NSum`;
  - d20: `FindRoot`, "non-real f".
- **Fewer digits at machine precision:**
  - c5, c6, c9: 2-D `NIntegrate` and minimax;
  - c7: sparse solve, 12 digits;
  - d6, d7: `NIntegrate`;
  - d9, d20: one digit short of 15.
- **Absent heads:**
  - d11: `FindRoot` on vector equations;
  - d4 `RiemannR`, d14 finite elements, d16 `WaveletPhi`, d17/d18 `NDEigen*`
    (each scored OK only because no other system answered either).
- **Slower:**
  - c3 SVD (2–9×);
  - c8 wp50 `NSum` (230×);
  - c10 machine `NSum` (7.7×);
  - d7 (3–5×);
  - d9 `MatrixPower` (4×);
  - d10 `FoldList` (9×);
  - c2 `NSolve` per bounce (2.6–2.9×);
  - d20 (3.9×).

### Final, v0.6460 (2026-11-09)

Full tables:
- `siam100_results/final_v0.6460.md`: experiments 166–225 (the originals and holdout A);
- `siam100_results/holdoutB_v0.6460.md`: holdout B.

Raw JSONs are in `raw_v0.6460.tar.xz`. Both runs were serial, on a quiet machine. WL was re-run for 171, 201 and 220 after their precision was aligned.

Originals plus holdout A, 676 cases:

| system | correct | inexact | wrong | none |
|---|---|---|---|---|
| **Mathilda 0.6460** | **484** | 127 | 26 | 39 |
| Python (scipy/mpmath) | 363 | 64 | 51 | 198 |
| wolframscript | 286 | 185 | 70 | 135 |
| Maxima | 214 | 25 | 6 | 431 |
| FriCAS | 127 | 14 | 14 | 521 |

On target:
- originals: 48/76 (v0.6377), then 64/76 (v0.6417), then **69/76**;
- holdout A: 262/600 (v0.6377, out of sample), then 438/600 (v0.6417), then **506/600**.

## Verdict

**Does the improvement generalise? Partly, and the improvement is substantial.** On the sealed holdout B (`SIAM100_HOLDOUT_B_RESULTS_v0.6460.md`), run once at v0.6460:
- Mathilda is on target in **395/600 (65.8 %)**. The out-of-sample figure before any fix was 262/600 (43.7 %, holdout A at v0.6377).
- Holdout A, now training data, stands at 506/600 (84.3 %).

So 133 of the 244 cases gained in sample carried over to unseen problems.

**Correctness.** On holdout B, Mathilda has the most correct answers of the five systems: 395, against 318 for Python, 275 for WL, 185 for Maxima and 108 for FriCAS. It gives 45 wrong answers, against 44 for Python and 61 for WL. 20 of those 45 come from one defect in D17 (every nodal-line length comes out `{0, 0}`).

**Speed.** Where Mathilda is correct, it is mostly also the fastest. The exceptions are:
- the D10 `FoldList` recurrences: 3–13× slower;
- C5 minimax: 10–20× slower;
- C3 SVD: 2–4× slower.

**The goal is not met.** The goal was "at least as good as every other system on every case". Mathilda meets it on two thirds of unseen cases.

**Next queue.** It is ranked by cases at stake, and the details are in the holdout-B results file. Fixing any of these makes holdout B training data, so the next measurement needs a holdout C.

| # | family | problem | cases |
|---|---|---|---|
| 1 | D17 | nodal length `{0, 0}`: one wrong-answer defect | 20 |
| 2 | D10 | `FoldList` random-recurrence speed | 20 + 14 in holdout A |
| 3 | C6, D8 | multi-dimensional `NIntegrate`: 50-digit 2-D/3-D integrals and infinite-range rational integrands time out | 18 + 7 in holdout A |
| 4 | D14 | FEM accepts only some linear operator forms (`femfail`/`femspec`) | 8 |
| 5 | D2 | `NSum` digit-length weights stop with `accgl` | 10 |
| 6 | D18, D19, D20 | unevaluated `NDEigenvalues` region forms, constrained design, functional roots | 18 |
| 7 | C3, C5 | SVD and minimax speed | — |
