# Known leak: per-call Solve/Reduce/FLINT/CAD leakage (surfaced via `Minimize`)

**Status:** open, pre-existing (NOT introduced by the v0.305–v0.307 `Minimize` work).
**Owner:** unassigned — pick up in a fresh session.
**Scope:** a `definitely lost` leak inside the symbolic `Solve`/`Reduce`/CAD/FLINT
machinery that grows **per call**, made visible by exact `Minimize` problems whose
optimum is an algebraic (irrational) number. It is the same class the v0.304
review called "pre-existing one-time Solve/Reduce/FLINT init leakage (~420
blocks)" — but it is **not** purely one-time; it scales with the number of heavy
algebraic `Reduce`/`Solve` calls.

## What was measured (valgrind `--leak-check=full`, macOS, this binary)

| Case | `definitely lost` |
|---|---|
| Startup only / rational-optimum `Minimize[{x^2+y^2,(x-2)^2+y^2==1},{x,y}]` (BASE) | 13,440 B / **420 blocks** |
| `Minimize[{x^2+y^2,(x-2)^2+(y-3)^2<=1},{x,y}]` ×1 (MIN, radical-cleared cert) | 20,112 B / **675 blocks** |
| same ×5 | 58,123 B / **1,907 blocks** |
| the raw `Reduce[... && u^2==13 && u>=0 && x^2+y^2 < 1/169(-39+3u)^2+1/169(-26+2u)^2, {x,y,u}, Reals]` ×1 | 16,072 B / **549 blocks** |

Per-call growth from the ×1 vs ×5 points: `(1907 − 675) / 4 ≈ 308 blocks / call`.
The BASE (rational-optimum) case leaked **nothing beyond the 420-block startup
baseline** — it stays rational, so it never exercises the leaky FLINT qqbar / CAD
paths. The leak appears only once the computation touches algebraic numbers
(`Sqrt[13]` via FLINT qqbar in candidate ordering, and the multi-variable CAD in
the lower-bound certificate).

## Diagnosis — it is NOT in `minimize.c`

The allocation stacks of the lost blocks pass through `minimize.c` **only as
call-sites** that invoke the evaluator:
- `mz_beval (eval.h:137)`  ← `eval_and_free` of a `Reduce`/`Solve` expression
- `mz_entails (minimize.c:368)`, `mz_run/mz_exact_poly (minimize.c:573)`

The **allocating** frames (the blocks actually never freed) are all in the
Solve/Reduce/CAD/FLINT machinery. Most frequent frames among the lost blocks:

```
786  eval.c:2739        636  eval.c:2219        526  solve (various)
379  expr.c:129         375  flint (various)    256  modular.c:94
223  flint_bridge.c:384 223  expand.c:352       204  internal.c:195
182  eval.c:1759        174  expr.c:194         136  match.c:1649
128  poly.c:1263        128  flint_bridge.c:1158 122 solvepoly.c:1684
116  solve.c:1070        79  simp_util.c:49
```

None of `minimize.c`'s own allocations (`mz_clear_radicals`, `mz_rad_walk`,
`mz_rad_push_*`, the `malloc`'d `vall`/`parts`/`cand` arrays, the `expr_copy`s in
`mz_prove_lower_bound`) own any lost block — they are all paired with frees.
So `Minimize` is a **witness**, not the cause: it issues more heavy algebraic
`Reduce`/`Solve` calls than most code, each of which leaks.

## How to reproduce

```bash
make -j
cat > /tmp/vg.m <<'EOF'
Do[Minimize[{x^2 + y^2, (x-2)^2 + (y-3)^2 <= 1}, {x, y}], {5}];
EOF
valgrind --leak-check=full --error-exitcode=0 ./Mathilda -file /tmp/vg.m 2>&1 \
  | grep 'definitely lost:'
# compare to a rational-optimum case that stays off the FLINT/CAD paths:
#   Minimize[{x^2 + y^2, (x-2)^2 + y^2 == 1}, {x, y}]   -> ~420 blocks (startup only)
```

To see the owning frames (confirm they are NOT in `minimize.c`):

```bash
valgrind --leak-check=full ./Mathilda -file /tmp/vg.m 2>&1 \
  | grep -E "definitely lost|flint|solve|eval\.c|expr\.c|minimize"
```

## Where to dig next (future session)

The per-call leak is a real bug in the shared symbolic engine, worth its own
campaign independent of `Minimize`. Suggested order, cheapest first:
1. **Isolate at the `Reduce`/`Solve` boundary.** Run valgrind on a bare
   `Reduce[...algebraic-coefficient system..., Reals]` and on `Solve[...,Reals]`
   that returns `Sqrt`/`Root` values, with `--num-callers=40`. The RAW-Reduce ×1
   row above (549 blocks) already shows one `Reduce` leaks ~a few hundred blocks.
2. **Start at the highest-frequency owning frames:** `eval.c:2739`, `eval.c:2219`
   (these look like evaluator-side temporaries not freed on some path), then the
   FLINT bridges `flint_bridge.c:384` / `:1158` and the CAD path
   (`solvepoly.c:1684`, `solve.c:1070`, `modular.c:94`). Check whether these are
   genuine leaks or `still reachable`/arena memory valgrind mislabels across the
   GMP/FLINT custom allocators (the `tc_gmp_alloc_*` machinery in `core.c`).
3. **Add a regression gate** once a fix lands: a valgrind step over a fixed
   `Reduce`/`Solve` script asserting `definitely lost` stays at the startup
   baseline across N iterations (the per-call delta must be 0).

## Why it is safe to ship the `Minimize` work now

`minimize.c`'s own memory is clean (verified above and by the `test_minimize.c`
memory-smoke loop). The `Minimize` feature does not add a new leak; it only
exercises more of a pre-existing engine-wide leak. Fixing that engine leak is a
separate, larger effort tracked by this file.
