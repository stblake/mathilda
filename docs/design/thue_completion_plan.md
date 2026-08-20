# Thue solver — completion plan (algorithmically clearing the declines)

**Purpose.** A self-contained roadmap to make `Solve[F(x,y)==m && Element[{x,y},
Integers], {x,y}, Integers]` solve the equations it currently declines, so a
fresh context can start executing without re-deriving anything. Read this first,
then `docs/references/thue/ALGORITHM_NOTES.md` (the implemented algorithm) and
`benchmarks/88-thue-equations/` (the verification gate).

---

## 1. Where we are

The Tzanakis–de Weger engine is **built, correct, and fast** for its scope:
monic irreducible binary forms `F(x,y)=m` with `|m|=1` over a **monogenic** real
field. Pipeline (all in place):

- `src/numbertheory/numberfield.{c,h}` (+`numberfield_internal.h`) — field setup,
  `arb`/`acb` embeddings, `disc`, **Gate 1** (Dedekind maximal-order cert),
  exact Sturm signature.
- `src/numbertheory/nfunits.{c,h}` — **Gate 2** fundamental units + regulator via
  small-norm search + **p-saturation**.
- `src/solvethue.{c,h}` — reduce to unit equations, Baker(Waldschmidt) bound,
  de Weger LLL reduction (`arb`, adaptive precision), Q-dependent case iii
  (relation-detection + L-trick), reconstruct + exact verify. `si_solve_thue`
  dispatches in `src/solveint.c`.

**Contract (must be preserved by everything below):** *provably complete or
DECLINE*. Never a wrong or partial answer. `benchmarks/88` (PARI/GP `thue()`
oracle over ~100 equations) is the gate: any `WRONG`/`CRASH` fails.

**Benchmark 88:** after M1 (Voronoi units) + M2 (general `m`) + M3 (Round-2
maximal order) + reducible forms (§6) + Minkowski-LLL O_K basis + M2b (rank-2
`|m| != 1`) + M5 (totally complex), **99 CORRECT / 5 DECLINE / 0 WRONG / 0 CRASH**
(48/56 before M1, 56/48 after M1, 65/39 after M2, 81/23 after M3, 94/10 after
reducible forms, 96/8 after Minkowski-LLL, 98/6 after M2b). Of the 5 declines, 2
are the reducible perfect powers PARI also refuses (correct on both sides); the 3
genuine gaps are all **M4** (rank ≥ 2 unit finding): the large-regulator
non-monogenic quartic `Q(10^{1/4})` (`±1`) and the quintic `x^5-5y^5`. Totally
complex (`r1=0`, M5) is now DONE. Cross-checked vs PARI over a 270-case
`|m| != 1` grid (M2), a 130-case non-monogenic grid (M3), a 90-case reducible
grid, a 150-case two-quadratic grid, and the checked-in randomized `grid.py`
(seed 20260820, 400 cases): 0 WRONG (1 case is a PARI `thue()` incompleteness,
adjudicated `PARI_WRONG`, Mathilda brute-verified correct).

### The 56 declines, by root cause (exact counts from the benchmark)

| # | bucket | why we decline | fix (this doc) |
|---|---|---|---|
| ~28 | **Non-monogenic field** (Gate 1) | we only certify `Z[θ]=O_K` | §3.2 Round-2 maximal order — ✅ DONE cubics+quartics (M3) |
| 11 | **Large-regulator units** (Gate 2) | a coeff-box search can't find units with large coordinates (ℚ(∛41): 24-digit coords, reg 56) | §3.1 Voronoi — ✅ DONE for complex cubics (M1) |
| 14 | **`\|m\| ≠ 1`** | reduction assumes `β=x−θy` is a *unit* | §3.3 μ-enumeration — ✅ DONE for rank-1 complex cubic (M2) |
| 1 | **Totally complex** (`r1=0`) | Siegel setup needs a real `i0` | §3.4 elementary \|Im\| bound — ✅ DONE (M5) |
| 2 | reducible perfect powers | *not a Thue equation* (infinitely many pts) | already correct (decline; PARI errors too) |

---

## 2. The shared prerequisite: grow the number-field layer toward `bnfinit`

Buckets 2, 3, and 1 are three views of the same missing capability — Mathilda's
number-field layer stops at "monogenic order + small units", whereas PARI's
`bnfinit` gives the **maximal order, the full unit group, and ideal
arithmetic**. Build these as reusable pieces in `src/numbertheory/`, then each
bucket is a thin application. FLINT/ANTIC gives exact `nf_elem` field
arithmetic but **none** of these three, so they are from-scratch (well-specified
classical algorithms; references below).

Three core pieces, in dependency order:

1. **Maximal order `O_K`** (integral basis as an HNF ℤ-basis `{ω_1..ω_n}`) —
   Round 2 / Pohst–Zassenhaus. Prereq for everything: units live in `O_K^×`, and
   a non-monogenic field hides units in `O_K ∖ Z[θ]`.
2. **Fundamental units + regulator** over `O_K` — Voronoi (cubic) / Buchmann
   (general). Replaces the coeff-box search's regulator ceiling.
3. **Ideal arithmetic** — prime-ideal factorisation of `(m)`, and a
   principal-ideal test with generator. Only bucket 1 strictly needs the *full*
   version; a lighter μ-enumeration (see §3.3) sidesteps most of it.

Everything stays exact / rigorously bounded (arb balls; every ACCEPT re-checked
in ℤ) to preserve the contract. **When any step can't be certified within
budget → DECLINE**, exactly as now.

---

## 3. The buckets

### 3.1 Large-regulator fundamental units (11 declines) — ✅ DONE (cubic), 2026-08-19

> **Status: M1 shipped for the complex-cubic case.**
> `src/numbertheory/nfvoronoi.c` (`nf_voronoi_unit_cubic11`) implements Voronoi's
> chain of relative minima for a rank-1 complex cubic (signature `(1,1)`), wired
> as the fallback in `nf_fundamental_units` after the box search fails. The
> neighbour search rescales the real Minkowski column by `1/U` to isotropise the
> anisotropic `|sigma2|<1` box, LLL-reduces (`lll_reduce_q` + identity block for
> the transform), and enumerates; `theta` is tracked exactly in mpz and the
> certifier (`eval_coords_mod`/`p_saturate`) is now mpz end-to-end (no int64
> cliff). Clears the cubic part of this bucket: `Q(cbrt 15/41/42/97)` now solve
> (benchmark 88: CORRECT 48→56, 0 WRONG/0 CRASH; regulators match PARI). Rank-2
> **totally-real** cubics (Buchmann's 2-D generalisation) remain — but the
> "Thomas family" cases in benchmark 88 are *reducible* (`x=-1` is always a root),
> i.e. correct declines, and the genuine small-reg totally-real cubics are already
> box-found, so there is no validatable rank-2 target here yet (deferred to M4).

**Root cause.** `nf_fundamental_units` (`nfunits.c`) enumerates algebraic
integers `Σ c_i θ^i` with `|c_i| ≤ box` and keeps `|N|=1`. A fundamental unit's
coordinates grow with the regulator, so `ℚ(∛6)` (coord 6) is found but `ℚ(∛15)`
(coord 30) and `ℚ(∛41)` (11-digit, reg 56) are not. No coefficient — or even
T2-norm — enumeration reaches reg-56 units; PARI doesn't enumerate, it walks.

**Fix — Voronoi's algorithm for the fundamental unit(s).** Voronoi walks the
"chain of relative minima" of the ring lattice; its length is **polynomial in
the regulator** (not exponential in coordinate size), so reg 56 is a few hundred
cheap steps.

- **Signature (1,1) complex cubic** (rank 1, e.g. ℚ(∛d)): classical Voronoi. The
  fundamental unit is the product of the successive minima around one period of
  the chain. Implement `nf_voronoi_unit_cubic_11()`.
- **Totally real cubic (3,0)** (rank 2): the 2-dimensional Voronoi generalization
  (Buchmann's algorithm I). Two independent fundamental units from the 2-D chain.
- **Quartic and higher / rank ≥ 3**: Buchmann's generalized Voronoi (II), or the
  subexponential (Buchmann–Lenstra) relation method. Larger; do after cubics.

**Wire-in.** Keep the current small-norm search as the FAST PATH (it already
solves small-regulator fields in ~1 ms); call Voronoi only when the search fails
to certify within the box (i.e., in the current "grow the box" loop, replace the
final give-up with a Voronoi attempt). Feed the resulting units into the SAME
p-saturation certifier (`p_saturate`) so the contract holds unchanged.

**Files.** `src/numbertheory/nfunits.c` (new `nf_voronoi_*`), possibly a new
`src/numbertheory/nfvoronoi.c`. Reuse `nf_norm_int`, the log-embeddings, and
`lll_reduce_q`.

**References.** Cohen, *A Course in Computational Algebraic Number Theory*
(CCANT) §5.7 (Voronoi, cubic); Buchmann, "A generalization of Voronoi's unit
algorithm I & II", *J. Number Theory* 20 (1986) 177–191, 192–209 (already in the
Tzanakis–de Weger reference list, [8]); Williams–Cormack–Seah for pure cubics.
The paper's **Appendix I** (Billevič) is the totally-real-quartic unit method it
actually used — implement that for the quartic case.

**Tests.** `binom-cubic-d15/d41/d42`, `adv-big-coef-cubic` (all monogenic per
PARI `nfinit`); assert Mathilda's set == PARI's. Regulators cross-check against
PARI `bnfinit(f).reg`.

---

### 3.2 Non-monogenic fields (~28 declines) — ✅ DONE (cubics + quartics), 2026-08-20

> **Status: M3 shipped — Round 2 maximal order + O_K-basis unit search.**
> `nf_round2_maximal_order` (`src/numbertheory/nfround2.c`, FLINT matrices)
> computes `O_K` by Pohst–Zassenhaus (p-radical Frobenius kernel, ring of
> multipliers via HNF, iterate), returned as an integral basis `(1/D)W` + index +
> `d_K`. `nf_field_create` runs it instead of declining; the struct carries
> `(W,D,index)` (monogenic = `I,1,1`, unchanged). The unit search (`nfunits.c`)
> walks the O_K lattice `L={Σ c_i W[i]}` testing `|N(v)|=Dⁿ` (v the θ-numerator),
> with `embed/D` and p-saturation's `D⁻¹ mod q`; the reconstruction is unchanged
> (θ-coords). Clears the non-monogenic cubics (`Q(cbrt 10/12/17/19/20)`, the
> Dedekind cubic) and quartics (`Q(d^{1/4})`, incl. index-16 `Q(12^{1/4})`).
> Benchmark 88: CORRECT 65→81; 130-case PARI grid 0 WRONG. Validated: index +
> `d_K` match PARI `nfinit`. Remaining: larger-regulator non-monogenic quartics
> (need O_K-Voronoi, the M4 unit finders) and general `m` over non-monogenic
> (gated off — M2×M3 follow-on).

**Root cause.** Gate 1 (`nf_is_maximal_order`, `numberfield.c`) only accepts
`Z[θ]=O_K`. Pure quartics `ℚ(d^{1/4})`, the simplest-cubic family, `ℚ(∛17)`
(index 3) etc. have index > 1, so we decline at `nf_field_create`.

**Fix — Round 2 (Pohst–Zassenhaus) maximal order.** When Dedekind reports
`Z[θ]` non-p-maximal at a prime `p` (`p²|disc`), enlarge:

1. p-radical `I_p = { x ∈ O : x^{p^k} ∈ pO }` = kernel of the Frobenius power map
   on `O/pO` (an 𝔽_p-nullspace; reuse `src/linalg/nullspace.*`).
2. ring of multipliers `O' = (I_p : I_p)` as an HNF ℤ-module (reuse
   `linalg_hnf`).
3. replace `O ← O'`, verify `disc(O')·[O':O]² = disc(O)` (exact index drop),
   repeat until p-maximal. Iterate over all `p` with `p²|disc`.

Output: `O_K` as an HNF integral basis `{ω_1..ω_n}` and `d_K`.

**Downstream refactor (the real work).** Everything currently assumes the
`θ`-power basis `Z[θ]`:
- **Unit search** must enumerate `O_K` elements `Σ c_i ω_i` (integer `c_i`), not
  `Z[θ]` — that's how it finds the units in `O_K ∖ Z[θ]`.
- **Reconstruction** in `solvethue.c`: `β = x − θy` still lies in `Z[θ] ⊆ O_K`,
  so its `O_K`-coordinates are fine; the "is it of the form `x−θy`" test converts
  `O_K`-coords → `θ`-basis and checks the `θ²..θ^{n-1}` coords vanish. Keep FLINT
  `nf_elem` over `ℚ(θ)` for arithmetic (it's basis-agnostic); use the `O_K` basis
  only where integrality / the unit lattice matters.
- **Signature/embeddings/Baker constants** are field invariants — unchanged.

Keep the monogenic fast path (`Z[θ]=O_K`) exactly as now; Round 2 only runs when
Dedekind fails.

**Files.** `src/numbertheory/numberfield.c` (`nf_round2_maximal_order`, return an
integral basis in the struct), `nfunits.c` (search over the `O_K` basis).

**References.** Cohen CCANT §6.1 (Round 2), Thm 6.1.3; Pohst–Zassenhaus,
*Algorithmic Algebraic Number Theory* Ch. 5.

**Tests.** `binom-quartic-d{2,3,5,7,10}` (`ℚ(d^{1/4})` non-monogenic),
`thomas-t*`, `nonmono-d17/d19`; set == PARI.

---

### 3.3 General `m` (`|m| ≠ 1`, 14 declines) — ✅ DONE (rank-1 complex cubic), 2026-08-19

> **Status: M2 shipped for the complex-cubic case (route a).**
> `thue_norm_reps_cubic11` (`src/solvethue.c`) enumerates the bounded-norm
> representatives μ: a canonical orbit rep reduces to `|σ1(μ)| ∈ [1, e^R)`, the
> norm forces `|σ2(μ)| ≤ √|m|`, and those embedding bounds map through the
> inverse Vandermonde to a finite coordinate box (kept if `N == m`; over-coverage
> is safe — `solbuf` dedups). `thue_exponent_bound` is now μ-aware: `δ` gains the
> `log|μ^(k)/μ^(j)|` term (looped over μ), and `C4`/`Y2p`/`V0` use the M-set
> min/max — all **over-estimates** (a too-large bound is safe; a too-small one
> would miss solutions). Cross-checked vs PARI over a 270-case grid: 0 WRONG.
> Clears `binom-cubic-d2-m{2,3,4,5,9,10,73,100}` + `nosol-d2-m4` (benchmark 88:
> CORRECT 56→65).
> **M2b (2026-08-20): rank-2 totally-real `|m|≠1`.** `thue_norm_reps_cubic11`
> now also enumerates the rank-2 case — a rep reduced into the fundamental
> parallelogram of `<L(ε1),L(ε2)>` has, per real embedding `i`, `|L_i(µ)| ≤
> |L_i(ε1)|+|L_i(ε2)|+|log|m||/3`, giving a finite coordinate box. Clears
> `cyclic-cubic-m2`, `nosol-cyclic`; `x³−3xy²+y³=8 → 6 pts`; 48-case rank-2 PARI
> grid 0 WRONG (benchmark 88: CORRECT 96→98). General `m` over a non-monogenic
> field still declines (M2×M3 follow-on).

**Root cause.** `thue_enumerate` bails at `if (mpz_cmpabs_ui(m,1)!=0) goto done;`.
The reduction assumes `β=x−θy` is a **unit**. For general `m`, `β` is an
algebraic integer of **bounded norm** `N(a0·x − θy) = a0^{n-1}·m`, so
`β = μ · ∏ εₖ^{bₖ}` where `μ` ranges over a **finite** set of bounded-norm
representatives (one per (ideal class × orbit)).

**Fix — enumerate the μ, then loop the existing reduction per μ.** This reuses
the whole current engine once the μ-set is in hand; the μ-set is the only new
piece. Two ways to get it (prefer the first — it avoids full ideal machinery):

- **(a) Direct bounded enumeration (Tzanakis–de Weger §II.1, recommended).**
  Once fundamental units are known (§3.1), a μ with `N(μ)=a0^{n-1}m` can be
  **size-reduced** by units so its log-embedding lies in the fundamental domain
  of `⟨log εₖ⟩`. That bounds `μ`'s coordinates → a finite Fincke–Pohst
  enumeration of `O_K` elements with `N=a0^{n-1}m` in that domain. Each survivor
  is a μ; dedupe unit-equivalents.
- **(b) Ideal-theoretic.** Factor `(a0^{n-1}m)` into prime ideals (from the
  splitting of the rational primes dividing it), enumerate integral ideals of
  that norm, test principality, extract a generator (needs the class group +
  `bnfisprincipal`-equivalent). More general but much heavier — defer.

Then: for each μ, run `thue_exponent_bound` / reduction / enumeration with the
Siegel form's `δ` including `log(μ^{(k)}/μ^{(j)})` (the code already has the
`μ`-ratio slot — `α_0` in `ALGORITHM_NOTES.md` — it just currently uses `μ∈{±1}`).
Union the per-μ solution sets; every candidate is still exact-verified.

**Files.** `src/solvethue.c` (drop the `|m|=1` gate; add `μ`-enumeration and the
per-μ loop; generalise the `δ` term), reusing §3.1 units.

**References.** Tzanakis–de Weger §II.1 (the `μ`, `M` set); Smart, *The
Algorithmic Resolution of Diophantine Equations*, Ch. IV; Cohen CCANT §4.8, §6.5
(ideals) for route (b).

**Tests.** `binom-cubic-d2-m{2,3,10,73}` (nonempty) and `-m{4,5,9,100}` (proved
`{}`), `cyclic-cubic-m2`; set == PARI.

---

### 3.4 Totally-complex fields (`r1 = 0`) — ✅ DONE (2026-08-20), by a simpler route

> **Status: M5 shipped — the elementary imaginary-part bound.**
> The plan below proposed torsion enumeration + a complex-`i0` Baker bound. None
> of it was needed: the totally-complex structure gives a **direct, elementary,
> rigorous** bound. Every root `theta_i` is non-real, so for real integers
> `x, y`, `|x − theta_i y| ≥ |Im(theta_i)|·|y|` — no factor can be small (there
> is no real root for `x/y` to approach), hence
> `|m| = |a0| ∏_i |x − theta_i y| ≥ |a0|·|y|^n·∏_i |Im theta_i|` and
> `|y| ≤ (|m| / (|a0|·∏|Im theta_i|))^{1/n}`. No units, no torsion, no Baker/LLL.
> `thue_solve_totally_complex` (`src/solvethue.c`, wired at the old `r1<1`
> decline) computes `∏|Im|` as a certified `arb` lower bound (→ a rigorous
> over-estimate of `|y|`) and closes each `y` by exact univariate root-finding.
> Solves the WHOLE family, **any `m`**: `Q(ζ5)` cyclotomic quartic (6 pts),
> `x^4+y^4 == {1,2,17,82}` and `== 3 → {}`, `Phi_7`/`Phi_10`. Benchmark 88 CORRECT
> 98→99 (`quartic-cyclotomic`, the last totally-complex decline). The grid also
> surfaced a PARI `thue()` incompleteness on another `Q(ζ5)` generator
> (`…==5`, PARI `[]`, true set `{(1,2),(−1,−2)}`) — now adjudicated (grid
> `PARI_WRONG`, Mathilda brute-verified correct).

**Original root cause (superseded).** `thue_exponent_bound` requires `s = r1 ≥ 1`
(a real embedding for the type index `i0`), and the enumeration uses torsion
`{±1}`. The torsion / complex-`i0` port became moot once the `|Im|` bound made the
Baker machinery redundant for `r1 = 0`.

**References.** Tzanakis–de Weger Lemma 1.1 (the `i0>s` branch) — for context; the
shipped method is the elementary geometry-of-numbers bound above.

---

## 4. Recommended order & milestones

Each milestone is independently shippable and gated by benchmark 88 (0
WRONG/CRASH, coverage strictly up).

1. **M1 — Voronoi units, cubic** (§3.1 cubic): ✅ **DONE (2026-08-19)** — complex
   cubic (rank 1) shipped in `nfvoronoi.c`; `Q(cbrt 15/41/42/97)` solve, bench 88
   CORRECT 48→56. Totally-real rank-2 (Buchmann 2-D) rolled into M4 (no target yet).
2. **M2 — General `m` via μ-enumeration** (§3.3, route a): ✅ **DONE (2026-08-19)**
   for rank-1 complex cubic — μ-enumeration + μ-aware Baker bound; bench 88
   CORRECT 56→65, 270-case grid 0 WRONG. Rank-2 totally-real `|m|≠1` deferred (M2b).
3. **M3 — Round-2 maximal order + `O_K`-basis unit search** (§3.2): ✅ **DONE
   (2026-08-20)** for cubics + quartics — `nfround2.c` + O_K-lattice unit search;
   bench 88 CORRECT 65→81, 130-case grid 0 WRONG.
4. **M4 — Voronoi/Billevič + Buchmann units, quartic & rank ≥ 3** (§3.1 rest):
   the remaining bucket-3 quartics/quintics (`Q(10^{1/4})`, `x^5-5y^5`), and
   robustness for M3's fields. **The only remaining coverage milestone.**
5. **M5 — Totally complex (`r1=0`)** (§3.4): ✅ **DONE (2026-08-20)** — the
   elementary `|Im|` bound, not the planned torsion/complex-`i0` port; solves the
   whole family for any `m`. Bench 88 CORRECT 98→99.

Rationale: front-load the tractable, high-coverage, low-refactor items (M1, M2)
before the big `O_K` refactor (M3). Units (M1/M4) are the shared foundation, so
they come first.

---

## 5. Verification (the gate for every milestone)

- **`benchmarks/88-thue-equations/`** is the oracle harness. After each
  milestone: `python3 run.py` must show **0 WRONG / 0 CRASH** and CORRECT strictly
  increased. `results.json`/`REPORT.md` track the arc. PARI/GP `thue()` and
  `bnfinit` are the cross-checks (monogenicity, units, regulators, solution
  sets). `python3 grid.py` is the **randomized** counterpart — a
  deterministic-seeded fuzz of random forms (degree 3–6, mixed `m`) vs PARI,
  which catches a wrong *finite* answer on forms nobody curated (seed
  `20260820`, 400 cases: 261 CORRECT / 0 WRONG). Run it too on any change that
  touches the solve/enumeration paths.
- Add the newly-covered families to `tests/test_thue.c` (rigorous auto-bound
  path) as pinned regressions.
- Per change: `make check-c99`, macOS `leaks` (MSL) = 0, all suites pass, and the
  `USE_FLINT=0` degrade build compiles.
- **Contract audit:** every new ACCEPT (a unit is fundamental, an ideal is
  principal, a μ is complete, an order is maximal) must be exact or a rigorous
  arb bound; where it can't be certified within budget, DECLINE. Grep the new
  code for any float-only decision that gates completeness.

---

## 6. Risks & notes

- **Scope creep toward re-implementing `bnfinit`.** Buckets 2/3/1 really are a
  maximal-order + unit-group + ideal layer. Keep each milestone minimal and
  benchmark-gated; resist building the general subsystem before a bucket needs
  it. (Round 2 and Voronoi-cubic are classical and bounded; the subexponential
  class-group method is the one to avoid until forced.)
- **The `O_K`-basis refactor (M3)** touches the unit search and the
  reconstruction test — do it behind the monogenic fast path so nothing already
  working regresses.
- **Performance already-earned** (benchmark 88 perf pass): power-precompute
  enumeration, unit-search certify-in-loop, p-saturation rank-stall early-out,
  adaptive reduction precision. Don't regress these; the new unit finders should
  be *fallbacks* after the fast small-norm search, not replacements.
- **Reducible forms** — ✅ **DONE (2026-08-20).** `thue_solve_reducible_form`
  (`src/solvethue.c`) factors `F(x,1)`; when it splits into ≥2 coprime factors
  the set is finite: enumerate the signed divisor assignments `∏ G_i(x,y)^{e_i} =
  m/content` and solve each system `{G_i=d_i}` (parametrise a linear factor's
  line, substitute into a second factor, integer-root the univariate, verify).
  Works for any `m`; 30-case PARI grid 0 WRONG; clears `thomas-t*`, `x^3-y^3=m`,
  `x^4-y^4=m`, biquadratic. A pure power of one factor (`(x-y)^3=1`, infinite /
  PARI-refused) DECLINEs. The **no-linear-factor** case (e.g. two irreducible
  quadratics `(x^2-3xy+y^2)(x^2-4xy+y^2)`) is handled by eliminating `y` via the
  resultant `Res_y(G_0-d_0, G_1-d_1)` (sampled + integer-interpolated), whose
  integer roots give the candidate `x`; 150-case 2-quadratic PARI grid 0 WRONG,
  clears `adv-close-roots-1`.
