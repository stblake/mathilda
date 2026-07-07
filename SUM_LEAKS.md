# Pre-existing memory leaks in the `Sum[]` cascade

Recorded 2026-07-02 while adding the infinite-`Product[]` families (which reuse
`Sum[]` internally). These leaks are **pre-existing** — they live in the older
`Sum` sub-algorithm stages, not in any of the code added for the product work.

## What leaks

Every infinite `Sum` that flows through the cascade leaks a handful of small
`Expr` nodes, even when it evaluates to the correct closed form. Measured with
macOS `leaks --atExit` (baseline `1+1` is **0 leaks**, so these are real):

| Input | Leaks / bytes |
|-------|---------------|
| `Sum[1/k^2, {k,1,Infinity}]`            | 14 / 800  |
| `Sum[1/(k(k+1)), {k,1,Infinity}]`       | 47 / 2592 |
| `Sum[Log[k]/k^2, {k,1,Infinity}]`       | 26 / 1424 |
| `Sum[Log[k]/k^3, {k,1,Infinity}]`       | 33 / 1808 |
| `Sum[1/k + Log[(k-1)/k], {k,2,Infinity}]` | 54 / 2896 |

The leak count grows with the structural size of the summand: harder inputs make
the earlier stages do more (and leak more) work before a later stage closes the
sum.

## Where it leaks

Every leaked block's allocation stack is rooted in the cascade dispatch, at the
`try_def(...)` calls for the stages that **run and bail before** the closing
stage:

```
builtin_sum (sum.c)
  -> dispatch_def  sum.c:169   (Sum`Gosper)
  -> dispatch_def  sum.c:171   (Sum`Trigonometric)
  -> dispatch_def  sum.c:172   (Sum`Rational)
     -> evaluate -> ... -> builtin_plus / builtin_times -> expr_new_function (calloc)
```

i.e. the intermediate `Expr` trees these stages build with `Together` / `Apart`
/ `Cancel` / `Solve` while probing an input are not fully freed on their bail /
fall-through paths. The same class is referenced by the memory note about a
"pre-existing GMP-rational leak" in the `NSum` numeric path.

## What is NOT leaking

Confirmed leak-clean (0 blocks trace into them):

- **All new `Product` stages** — `Product`Viete`, `Product`Cantor`,
  `Product`EulerPrime`, `Product`LogSum`, `Product`RationalInfinite`, and the
  infinite `Product`Geometric` extension. Every pure-product input above
  (`Product[(k^2-1)/(k^2+1),...]`, `Product[Cos[Pi/2^(k+1)],...]`,
  `Product[1+(1/3)^(2^k),...]`, `Product[1/(1-1/Prime[k]^2),...]`) reports
  **0 leaks / 0 bytes**.
- **The new `Sum` stages** — `Sum`LogZeta` and `Sum`LogRational`. No leaked
  block's stack passes through `sum_logzeta.c` or `sum_logrational.c`; the
  leaks on their inputs come entirely from the pre-existing earlier stages
  (Gosper / Trigonometric / Rational) that run first and fall through.

## Reproduce

```bash
make -j$(nproc)
printf '{"id":1,"expr":"Sum[1/k^2,{k,1,Infinity}]"}\n{"type":"quit"}\n' \
  | MallocStackLogging=1 leaks --atExit --list -- ./Mathilda
```

Group the `Call stack:` lines by the `dispatch_def sum.c:<line>` frame to see
which stage each block came from.

## Scope

Fixing these means auditing the bail/fall-through paths of `Sum`Gosper`,
`Sum`Trigonometric`, and `Sum`Rational` for owned `Expr*` trees (produced by
`sum_eval("Together"/"Apart"/"Cancel"/"Solve", ...)`) that are dropped without
`expr_free`. That is a self-contained cleanup of the existing `Sum` subsystem,
independent of the product work, and is left as a follow-up.
