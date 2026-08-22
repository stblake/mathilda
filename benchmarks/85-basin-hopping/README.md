# 85 — Basin Hopping (`NMinimize`, `Method -> {"BasinHopping", …}`)

Races Mathilda's `nm_basin_hopping` (`src/numerical_calculus/nm_basin_hopping.c`)
against [`scipy.optimize.basinhopping`][bh] on eight standard bounded multimodal
benchmarks. Both implement Monte-Carlo minimization (Wales & Doye 1997): each hop
is a uniform random displacement of the current point, followed by a local
minimization (the "quench"), accepted by a Metropolis rule on the two
locally-minimized energies, with an adaptive step size that targets a fixed
acceptance rate.

## Method

- **Same algorithm parameters on both sides** — scipy's defaults (`T 1`,
  `stepsize 0.5`, `interval 50`, `target_accept_rate 0.5`, `stepwise_factor 0.9`,
  `niter 100`), a **single run** (no multi-start) from a seeded random start in the
  box, and a fixed seed (`RandomSeed -> 1` / `seed=1`) so both runs are
  reproducible. The RNG streams differ (Mathilda SplitMix64, scipy PCG64), so the
  trajectories differ; the race is like-for-like on the algorithm, not the stream.
- **Check** = the objective at the global optimum, `Round[10^6 · f]` on the
  Mathilda side and `int(floor(1e6·f + 0.5))` on the scipy side. All eight cases
  agree to 1e-6 (0 CHECK-FAIL).
- 2-D objectives are written with bare variables so they hit Mathilda's compiled
  machine-precision path — the fair analogue of scipy's numpy objective. The 5-D
  Ackley is built over fresh global symbols and passed through `Evaluate[]` (the
  same idiom as `79-shgo` / `81-dual-annealing`). An inline `Sum[…x[i]…]` would fall
  to the interpreter and is ~1000× slower here — a Mathilda indexed-variable
  limitation unrelated to the Basin Hopping engine, not a property worth timing.

## Cases

Himmelblau, Booth, Beale, Rosenbrock, six-hump camel, Ackley (2-D and 5-D),
sphere — eight functions whose global a **single** run of each system reaches and
polishes to agreement.

## Result

**8/8 AHEAD, 0 CHECK-FAIL.** Mathilda reaches the identical global on every case
and is **~100×–650× faster per solve**:

| case | Mathilda | scipy | speedup |
|------|----------|-------|---------|
| 07 Sphere 2D       | 0.103 ms | 66.9 ms | ~650× |
| 04 Rosenbrock 2D   | 0.588 ms | 302.3 ms | ~514× |
| 01 Himmelblau 2D   | 0.328 ms | 140.3 ms | ~428× |
| 03 Beale 2D        | 0.651 ms | 236.2 ms | ~363× |
| 02 Booth 2D        | 0.269 ms | 87.0 ms  | ~323× |
| 06 Ackley 2D       | 3.1 ms   | 1.00 s   | ~323× |
| 05 Six-hump camel  | 0.670 ms | 157.5 ms | ~235× |
| 08 Ackley 5D       | 15.6 ms  | 1.60 s   | ~103× |

The gap is per-evaluation cost: a basin-hopping run spends its time inside ~100
local minimizations (both systems do the same number), and Mathilda's objective is
compiled to bytecode where scipy calls back into Python for every evaluation the
quench needs.

## Fair-comparison envelope

These eight functions are ones where a single run of both implementations reaches
the same global, so the 1e-6 check is unambiguous. Basin hopping's effectiveness on
harder landscapes depends on the **quench**, and that is the one deliberate
difference between the two: scipy quenches with L-BFGS-B, whose aggressive first
step sometimes *overshoots across a basin boundary*; Mathilda quenches with its own
well-behaved BFGS / augmented-Lagrangian / integer-descent polish, which does not.
Three regimes fall outside the clean race:

- **Mathilda reaches the global where scipy's single run stalls.** On Rastrigin-2D
  at `seed=1`, the adaptive-step walk reaches `0` while scipy's single run stalls in
  a ring minimum near `0.995` (its L-BFGS-B quench does not, on that stream, escape
  the ring). Included as a CHECK-FAIL if raced, so it is documented here rather than
  timed.
- **scipy's overshoot escapes a local basin a single Mathilda run does not.** On the
  1-D quartic `x⁴ − 3x² − x` (basins at x≈1.30 global, x≈−1.13 local) and on
  Styblinski-Tang (four separated basins), scipy's single run reaches the global on
  most seeds, whereas a single Mathilda run crosses those wide gaps only through the
  0.5 displacement and so depends on the starting basin. `MaxIterations` does not
  help (the conservative walk is genuinely stuck); `"SearchPoints" -> K` multi-start
  does (quartic 7/8 seeds at K=4, Styblinski 4/8 at K=6 — and both nail it at the
  seeds the unit tests pin).
- **Constrained / integer, which scipy cannot express.** `scipy.optimize.basinhopping`
  is box-only; Mathilda's quench handles inequality, equality, disjunctive, and
  `Element[·, Integers]` constraints (e.g. `min x + y` on the disk `x² + y² ≤ 9` →
  `−4.2426`, and `min (x−3)² + (y+2)²` over integers → `0` at `(3, −2)`).

## Run

```
HPC_PYTHON=/usr/local/bin/python3.11 python3 benchmarks/run_all.py \
    --only 85 --system mathilda,python
```

(writes `REPORT.partial.md` / `results/<date>-partial.json` only — the canonical
weekly artifacts are left untouched; do not commit the `.partial` files.)

[bh]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html
