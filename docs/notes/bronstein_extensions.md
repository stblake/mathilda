# Bronstein extensions: the theory behind the completion campaign

This note holds every mathematical statement the Bronstein completion campaign
(`src/internal/bronstein/BRONSTEIN_PLAN.md`, "2026-11-14 -- Completion campaign")
relies on beyond a direct transcription of the sources:

- [B] M. Bronstein, *Integration of elementary functions*, JSC 9 (1990) 117–173;
- [B91] M. Bronstein, *The Risch differential equation on an algebraic curve*,
  ISSAC '91;
- [SI] M. Bronstein, *Symbolic Integration I*, 2nd ed., Springer 2005;
- [T] B. Trager, *Integration of algebraic functions*, PhD thesis, MIT 1984.

Part A gives derivations of results the sources contain, written out because the
implementation depends on their exact form (the off-by-one in A2 was a bug
before it was a derivation). Part B states results the sources do not contain,
with proofs. Every proof in Part B is reviewed adversarially before its phase
starts, and every claim carries the regression test that pins it.

Notation: `k` is a differential field of characteristic 0 with constant field
`C`, `θ` a monomial over `k` (primitive, `Dθ = η ∈ k`, or hyperexponential,
`Dθ/θ = η ∈ k`), `K = k(θ)`. For a polynomial `p`, `lc(p)` is its leading
coefficient and `deg p` its degree in θ.

---

## Part A. Derivations of known results, as implemented

### A1. Parametric logarithmic derivatives in a tower (`rc_param_log_deriv`)

**Problem.** Given `f, η ∈ k = Q(x)(t_0..t_{m-1})`, a Log/Exp tower with
algebraically independent monomials, and θ = exp(∫η) transcendental over `k`:
decide whether `f = m η + Dz/z` with `m ∈ Z`, `z ∈ k*`. (With η omitted: `m = 0`.)

**Claim.** Let `P` be the irreducible factors over Q of `den(f) den(η)`, as
polynomials in `x, t_0..t_{m-1}`, excluding the exponential monomials `t_j`. Let
`E` be the exponential monomials of `k`. Then `f = m η + Dz/z` has a solution iff
`f` lies in the Q-span of

    G = {η} ∪ {Dp/p : p ∈ P} ∪ {Dt_j/t_j : t_j ∈ E}

with every coefficient an integer. The coefficients are then unique, and
`z = c Π p^{e_p} Π t_j^{s_j}`.

**Proof.** Write `z = c Π p^{e_p} Π_{j∈E} t_j^{s_j}`, with `p` running over the
irreducibles of `Q[x, t]` other than the exponential monomials. Unique
factorisation in `Q[x, t]` permits this, because the monomials are algebraically
independent; a log monomial `t_l` is an ordinary irreducible here. Then
`Dz/z = Σ e_p Dp/p + Σ s_j Dt_j/t_j`.

- **Every `p` is normal at its top variable.** Take the highest variable `v` of
  `p`. If `v` is `x` or a log monomial, `Dv ∈ k_{<v}`. If `v` is an exponential
  monomial other than `p` itself, the only special irreducible is `v`. Either way
  `p ∤ Dp` in `k_{<v}[v]`.
- **So `Dp/p` has a simple pole along `p` with residue 1.** No other `Dq/q` has a
  pole along `p`. Those with the same top variable are coprime to `p`. Those with
  a lower top variable are constants for the partial-fraction decomposition in
  `v`. `Dt_j/t_j ∈ k_{<t_j}` has no pole along any `p`.
- **Hence `e_p ≠ 0` puts `p` in `den(m η + Dz/z) = den(f)`,** unless the pole
  cancels against `m η`, in which case `p | den(η)`. So every `p` with
  `e_p ≠ 0` lies in `P`, and the representation lies in the span of `G` with
  integer coefficients.
- **`G` is Q-linearly independent.** A relation, cleared to integers, says
  `D(θ^a Π p^b Π t_j^c) = 0`, so `θ^a Π p^b Π t_j^c` is a constant. For `a ≠ 0`
  this makes θ algebraic over `k`, contradicting transcendence. For `a = 0` it
  contradicts unique factorisation in `Q[x, t]`, the variables being
  algebraically independent.
- **So the span coefficients are unique.** The equation is solvable iff they
  exist and are all integers.

The independence of the tower is a hypothesis. A tower in which `E^(x + x^2)`
sits beside `E^x` and `E^(x^2)` breaks it. Since v0.6546 the tower builder
detects that (`RtTower.dependent`), and certificates on such a tower are demoted.

### A2. RdeBoundDegree for a primitive monomial (`rde_bound_degree_prim`)

**Setting.** We solve `a Dq + b q = c` for `q ∈ k[θ]` with `Dθ = η ∈ k`. Write
`q = Σ_{i≤n} q_i θ^i`, so `D(q_i θ^i) = Dq_i θ^i + i q_i η θ^{i-1}`. Let
`da, db, dc` be the degrees of `a, b, c`.

- **`db > da`.** The top term of `b q` is `lc(b) q_n θ^{db+n}`, and nothing else
  reaches that degree. So `n ≤ dc − db`.
- **`db = da`, with `α = −lc(b)/lc(a)` not of the form `Dz/z` (`z ∈ k*`).** The
  `θ^{da+n}` coefficient is `lc(a) Dq_n + lc(b) q_n`. It vanishes iff
  `Dq_n/q_n = α`. So `n ≤ dc − db`.
- **`db = da`, with `α = Dz/z`.** Now `q_n = c z` (c a constant) can cancel the
  top coefficient.
  - Then `deg(lhs) ≤ da + n − 1`, which may equal `dc`, so `n = dc − da + 1` needs
    nothing more.
  - A larger `n` also needs the `θ^{da+n−1}` coefficient to vanish. That
    coefficient is
    `lc(a)(Dq_{n−1} + n q_n η) + a_{da−1} Dq_n + lc(b) q_{n−1} + b_{db−1} q_n`.
  - With `q_{n−1} = z u`, and using `lc(a) Dz + lc(b) z = 0`, it becomes
    `lc(a) z Du + c (n lc(a) η z + a_{da−1} Dz + b_{db−1} z)`. This vanishes iff
    `−β = D(u/c) + n η`, where `β = (a_{da−1} Dz + b_{db−1} z)/(z lc(a))`.
  - That is a limited-integration problem with a unique constant `s` (two values
    would put θ in `k + C`).
  - So **`n ≤ max(dc − da + 1, s)`**, with `s` counted only when it is a positive
    integer.
  - The first implementation used `max(dc − da, s)`. It gave the false
    `NoSolution` for `y' + y/x = 2 Log[x]/x^2`, whose solution `Log[x]^2/x` has
    degree `dc − da + 1 = 2`. The regression test is in `test_bronstein_rde.c`.
- **`db = da − 1`.** The `θ^{da+n}` coefficient `lc(a) Dq_n` vanishes when `q_n`
  is a constant. The `θ^{da+n−1}` coefficient is then
  `lc(a)(Dq_{n−1} + n q_n η) + lc(b) q_n`, which vanishes iff
  `α = −lc(b)/lc(a) = D(q_{n−1}/q_n) + n η`. So
  `n ≤ max(dc − da + 1, s)` with `s = LimitedIntegrate(α, η)`.
- **`db < da − 1`.** If the top two coefficients both vanished, θ would lie in
  `k + C`. So `n ≤ dc − da + 1`.

`LimitedIntegrate(α, η)` against the one primitive θ asks whether `α` has an
antiderivative in `K = k(θ)`. Any such antiderivative is `z + sθ + const`, by
the degree argument above. It is decided by the field integrator
(`rc_limited_integrate_prim`). An answer carrying a new logarithm is still
reported as undecided; the ParamRischDE work of Phase 2 removes that.

### A3. RdeBoundDegree for a hyperexponential monomial (`rde_bound_degree_exp`)

Here `D(q_n θ^n) = (Dq_n + n η q_n) θ^n`, so degrees are preserved.

- **`da ≠ db`.** `n ≤ dc − max(da, db)`.
- **`da = db`.** The top coefficient `lc(a)(Dq_n + n η q_n) + lc(b) q_n` vanishes
  iff `−lc(b)/lc(a) = Dq_n/q_n + n η`. That is a parametric logarithmic
  derivative (A1), and `n = m` is unique. So `n ≤ max(dc − da, m)`.

### A4. RdeSpecialDenomExp (`rde_special_denom_exp`)

A τ-adic pole of order `−n` in `q` meets `a Dq + b q`. Its lowest coefficient is
`(a(0)(Dq_{−n} − n η q_{−n}) + b(0) q_{−n}) τ^{−n}`. When `ν_τ(b) = 0` this can
vanish, iff `−b(0)/a(0) = Dq_{−n}/q_{−n} − n η`: a parametric logarithmic
derivative with `m = −n`. Hence `n = min(ν(c) − min(0, ν(b)), m)` when A1 finds
`m`.

### A5. PolyRischDECancel{Prim,Exp} (`rde_polyrischde_cancel`)

Here `b ∈ k*` and `D q + b q = c`.

- **The homogeneous equation is solvable in `K`** iff `−b = Dw/w` with
  `w = z θ^m`, which by A1 means `−b = m η + Dz/z` (`m = 0` for a primitive θ).
  - If it is, `D q + b q = w D(q/w)`, so `q = w ∫ c/w + λ w`.
  - `λ` is fixed by the one coefficient that `q` must not have: a τ-pole when
    `m < 0`, or the degree when it exceeds `n`.
- **If it is not, the top-down pass is forced.** Each equation
  `Ds + (b + iη)s = lc(c)` in `k` has at most one solution, since a second would
  give a homogeneous solution. So a failed pass is a proof.
- **If the test is undecided,** a failed pass proves nothing, and it is marked.

### A6. WeakNormalizer over a tower (`rc_weak_normalizer`)

The source is [SI] §6.1, transcribed for any monomial.

- `f = A/(E d_1)`, with `d_1` the normal simple poles.
- The residues are the roots of `R(z) = Res_θ(d_1, A − z E Dd_1)`.
- `q = Π_n gcd(A − n E Dd_1, d_1)^n` over the positive integer roots `n`.
- The integer roots of `R` are common roots of all its coefficients as a
  polynomial in the tower variables. So they are the roots of the gcd of those
  coefficients over `Q[z]`, read off its linear factors over Q.

---

## Part B. Extensions beyond the sources (each proved before it is coded)

| | Result | Needed for | Status |
|---|---|---|---|
| E1 | The Risch DE `z' + f z = g` on a simple-radical curve `Y^m = G` over a transcendental tower `k(θ)`: order bounds at the places over `k[θ]`, at θ = ∞ and at the special places; reduction of the coefficient system to parametric Risch DEs in `k` | Phase 5 (radical anywhere; `Tower:radical-mid-tower`, `RDEAlgebraic:undecided`) | open |
| E2 | Transcendental Risch over an algebraic base `k_j(y)(t_{j+1}..)`: correctness of the recursion through the differential-field record, and Liouville's theorem for a mid-tower radical | Phase 5 | open |
| E3 | Torsion by specialization with the special places included: goodness conditions and injectivity in characteristic 0 | Phase 3 (`LogPart:special-genus`) | open |
| E4 | Symbolic parameters: the verdicts are generic (valid off a proper Zariski-closed set), and which zero tests may be generic | Phase 4 | open |
| E5 | Ideals of the `d = gcd(e, m)` ramified places above one special point | Phase 3 (`Divisors:multi-ramified`) | open |
| E6 | Any gap found during implementation | -- | -- |

Notes collected so far:

- **[B] pp.155–157 (η' a combination of logarithmic derivatives).** This needs
  `η' = Σ e_i γ_i'/γ_i` with constants `e_i`. With rational `e_i`, θ would be
  algebraic over `k`, so over the constant field Q the case cannot occur in a
  valid tower. It occurs for irrational or transcendental exponents
  (`x^Sqrt[2]`, `x^a`), so it belongs to Phase 4. Checked 2026-11-14 against the
  paper. The suspected false certificate in `brIntegrateInKLogLike` was refuted
  by the same reading: Def 2.21 allows logarithms of any element of `k`.
