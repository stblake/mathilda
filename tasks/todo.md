# Diophantine: Solve::svars warning (A) + Booker cube-root-mod-d (B)

Context: `Solve[x^3+y^3==d+z^3 && ... , {x,y,z,y}, Integers]` stayed unevaluated
because (1) the var list `{x,y,z,y}` drops `d` (a free symbol → classifier
declines) and (2) even corrected, the box exceeds the enumeration/MITM budget.
The engine already implements the classical divisor method
(`si_two_power_solve`), which is fast for FIXED-target forms (`x^3+y^3==1729` in
72 ms at bound 1e5). Two agreed improvements:

## Part A — `Solve::svars` diagnostic  (low risk)
Goal: when a symbol appears in the system (esp. with its own constraints) but is
NOT among the solve variables, emit a warning instead of silently declining.

- [ ] In `solveint_solve_integer` (src/solveint.c), after parsing the var list,
      scan flattened conjuncts for "free" symbols: EXPR_SYMBOL leaves that are
      (a) not in the solve-var list, (b) not a Protected constant / known head
      (Pi, E, True, False, Integers, ...), (c) actually a bound-eligible atom
      (appear in an Equal/inequality position, not a function head).
- [ ] If any found, `fprintf(stderr, "Solve::svars: Equations may not give "
      "solutions for all \"solve\" variables.\n")` ONCE (dedupe). Matches the
      existing stderr message style (Clip::ncompl, Power::infy).
- [ ] Emit the warning but still proceed (do not change the return value); the
      warning fires on the path that currently returns NULL for this shape.
- [ ] Keep it low-false-positive: only warn when the free symbol appears inside
      an (in)equality constraint (strong signal it was meant as a variable).
- [ ] Verify: the user's exact line now prints the warning; legitimate
      parametric solves (`Solve[a x==b,x]`) do NOT spuriously warn.

## Part B — Booker-style cube-root-mod-d for `x^3+y^3+z^3 == k`  (higher risk)
Goal: lift reach for the FIXED-k three-cubes problem from ~1e6 coords toward
~1e7–1e8 interactively (10–100x), by replacing "enumerate z, factor k−z^3" with
"enumerate d, cube-root k mod d". NOT aiming for Booker's 1e16 (that needs the
factorless sieve + cluster). Soundness is paramount: exact verify every hit;
DECLINE (stay unevaluated) rather than return an incomplete set.

Math (from Booker, "Cracking the problem with 33"):
  k − z^3 = x^3 + y^3 = (x+y)(x^2−xy+y^2);  d=|x+y| divides |k−z^3|;
  s = sgn(k−z^3)·d;  disc = (4|k−z^3|/d − d^2)/3;  x,y = (s ± sqrt(disc))/2.
  d | (k−z^3)  ⟺  z^3 ≡ k (mod d).

- [ ] Trigger detection: exactly 3 variables, each appearing only as v^3 with
      coefficient +1, constant term = −k, k a nonzero integer; fully box-bounded.
      (Only all-+1 so all three pairings are symmetric — keeps completeness
      reasoning clean. Other sign patterns fall back to existing method.)
- [ ] New helper `all_cube_roots_mod(k, d, roots_out)`: ALL z in [0,d) with
      z^3≡k (mod d). Factor d (df_factor_mpz); per prime power find all roots
      (p≡2 mod 3: unique via inverse-of-3 exponent; p≡1 mod 3: one root ×
      {1,ω,ω^2} where ω from sqrt(−3); p=3 and gcd(k,d)>1 handled explicitly);
      CRT product of all combinations. Correctness is the crux — mirror the
      proven logic in powermod.c (rth_root_mod_pe) but return ALL roots.
- [ ] Enumerate d in [1, D_max], D_max = min(alpha*B, budget). For each root
      class, walk z in the AP within the window, compute disc, test perfect
      square, recover x,y, bounds-check, `si_verify` exactly, emit.
- [ ] Role-loop over which variable is the "modular" one (3 passes) so the
      "two-largest" pairing (Booker's d<alpha·B bound) is always covered;
      dedupe solutions (canonical sorted tuple).
- [ ] Degenerate families: k−z^3==0 (z=cbrt(k), y=−x) and the y=z Thue slice —
      catch via a cheap direct small-coordinate check so nothing is missed.
- [ ] Completeness/decline: only RETURN a list when the box is provably covered
      by D_max (alpha-bound); else DECLINE. Budget-gate like the existing paths.
- [ ] Wire into the special-forms dispatch BEFORE the enumeration size guard.

### Verification (B is only shippable if these pass)
- [ ] Cross-check vs the EXISTING divisor method on overlapping boxes: identical
      solution SETS for many random k over a small box (independent code paths).
- [ ] Known solutions: 1729-style, x^3+y^3+z^3==29 → {1,1,3}; a constructed k
      with a ~1e7 coordinate found where the old path declines.
- [ ] `make check-c99`; clean `gcc -std=c99 -Wall -Wextra`; valgrind clean on a
      representative run; no packed-array surfaces affected (symbolic solver).

## Review

**Part A (`Solve::svars`) — DONE.** `si_warn_free_symbols` in `src/solveint.c`
scans only inequality/ordering conjuncts for symbol atoms not in the solve-var
list and not `Protected` (so operator heads + constants like Pi/E are skipped),
emits `Solve::svars` once (deduped by `expr_hash` of the system against the
evaluator's fixed-point confirm re-entry). Verified: the user's `{x,y,z,y}` query
now warns once; all-vars-present and bare-parameter-in-equation cases stay silent.

**Part B (Booker three-cubes) — DONE.**
- `si_all_cube_roots_mod` (+ `si_croots_mod_p` Pohlig-Hellman in the 3-Sylow for
  p≡1 mod 3, brute for prime powers, CRT product). Exposed as
  `Solve\`CubeRootsMod[k,d]`; **verified against brute force on 94 430 (k,d)
  pairs, 0 mismatches** (incl. p≡1 mod 9, p=3^e, p|k, composites, large primes).
- `si_solve_three_cubes_booker`: part1 (small-coord + divisor-solve, covers the
  `(a,-a,∛k)` family) ∪ part2 (Booker α-bound over 3 roles), verify-every-hit,
  decline when a box would exceed `SI_BK_SOLCAP`. **Cross-checked vs exhaustive
  box enumeration over 801 targets k, 0 mismatches.**
- Reach demo (no force): `x^3+y^3+z^3==2` over `[-200000,200000]^3` (2B=4e5, which
  the classical path declines) → 195 complete/verified solutions incl.
  `(162001,-161999,-5400)` in **0.39 s**.
- Bug found+fixed along the way: `build_result` used an O(n²) selection sort that
  hung on large solution families → replaced with O(n log n) `qsort`.

**Verification:** `solve_integers_tests`, `solve_tests`, `solve_corpus_tests`
(0/97 non-PASS) all green; `make check-c99` PASS; no valgrind leaks attributable
to new code (only macOS libobjc/dyld baseline noise). Docs: changelog
`2026-08-17.md` + `solutions-of-equations.md` updated.

**SPF sieve + 128-bit reach — DONE (follow-on).** `si_build_spf` (per-solve
smallest-prime-factor table, O(log d) factoring) + `__int128` part-2 arithmetic
(`si_isqrt_i128`) + `SI_BK_DMAX` 3e5→3e6 + `SI_BK_MAXNODES=1e9` candidate
backstop. Reach ~1e6 → ~1e7 coords: `x^3+y^3+z^3==2` finds the point at
**5 821 795** (radius-6e6 box) in ~16 s (radius 2e6 → 1 971 055 in ~4 s), all
verified; oversized boxes decline instantly (no hang). Correctness re-checked:
801-target box cross-check and 45 149-pair primitive check both 0 mismatches;
`make check-c99` clean; `solve_integers_tests` pass.

**Known scope limits (deliberate, documented):** Booker gated to |k|<~1e9 and
Dmax≤3e6 (coords ~1e7). Perfect-cube / large-family boxes decline rather than
emit O(B) tuples. Not Booker's 1e16 (needs batch inversion + a cluster).
`MATHILDA_BK_FORCE=1` bypasses the size gate for validation.
