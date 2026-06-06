# Known Memory Leaks

A running log of **pre-existing** memory leaks discovered in the codebase but
not yet fixed (because the fix is out of scope for the change that surfaced
them). Each entry should give enough detail to fix it later without re-doing
the investigation.

How to reproduce a leak check on macOS (valgrind reports a constant ~12,800 B /
400 blocks of dyld/libobjc/Accelerate init noise — always diff against a
`Sin[1.0]` baseline and only count stacks that reference `src/` files):

```bash
printf 'Sin[1.0]\nQuit[]\n'      | valgrind --leak-check=full ./Mathilda 2>baseline.txt
printf '<expr>\nQuit[]\n'        | valgrind --leak-check=full ./Mathilda 2>run.txt
grep "definitely lost:" baseline.txt run.txt        # compare totals
grep -E "\.c:[0-9]+\)" run.txt                       # inspect src frames
```

---

## `builtin_inequality`: residual `Inequality` arg buffer leaked

- **Where:** `src/comparisons.c:436`, in `builtin_inequality`.
- **What:** when an inequality chain has *undecided* pairs (e.g. `0 < p < 1`,
  `3 < a < 5` with symbolic variables), the function builds a residual
  `Inequality[...]` from a `malloc`'d `out_args` buffer (allocated at
  `comparisons.c:408`). `expr_new_function` **copies** that array into its own
  storage, so the caller must `free(out_args)` afterward. The single-surviving-
  pair path (`out_n == 3`, line 433) does free it, but the general
  residual-Inequality return at line 436 does **not** — leaking the buffer
  (`cap * sizeof(Expr*)`, i.e. 32 bytes for a simple chain, scaling with the
  number of undecided pairs).
- **Trigger:** evaluating any partially-undecided inequality. Surfaced by
  `PowerExpand[..., Assumptions -> 0 < p < 1]` etc., but the leak is in the
  evaluator's handling of the inequality itself, not in `PowerExpand` —
  confirmed by reproducing it with the bare expressions `0<p<1`,
  `3<a<5 && -2<b<-1 && 7<c<9`, `z<0` and no `PowerExpand` (identical
  +160 B / +4 blocks vs the `Sin[1.0]` baseline).
- **Stack (representative):**
  ```
  malloc (...)
  builtin_inequality (src/comparisons.c:408)
  evaluate_step (src/eval.c:812)
  evaluate (src/eval.c:1079)
  ...
  ```
- **Likely fix:** add `free(out_args);` before `return` at line 436 (the array
  elements are owned by the new `Expr`; only the array buffer needs freeing —
  same idiom as the `out_n == 3` branch and as `src/expand.c`). A quick audit of
  the other early-return paths in this function for the same pattern is
  warranted.
- **Severity:** low (small, bounded per call; not in a hot loop), but it scales
  with how often symbolic inequalities are evaluated.
- **Discovered:** 2026-06-06, during `PowerExpand` implementation.

---

## `PossibleZeroQ`: Schwartz–Zippel numeric stage leaks on symbolic radicals

- **Where:** `src/zero_test.c`, in the Stage-3 Schwartz–Zippel path
  `decide_schwartz_zippel` (zero_test.c:660) → `decide_numeric` (556) →
  `evaluate_rung` (529) → `numericalize_at` (→ `numericalize` in `src/numeric.c`).
  In an `-O3` build the whole chain is inlined, so `leaks` roots every block at
  `zero_test_decide` (zero_test.c:720, the `return decide_schwartz_zippel(e);`)
  with the allocation in `expr_new_function` (expr.c:81).
- **What:** when `PossibleZeroQ[e]` is given a symbolic expression that the
  structural (`decide_structural`) and rational (`decide_rational`) stages can't
  settle — characteristically anything containing a **radical of a symbol**
  (`x^(1/3)`, `Sqrt[x+1]`, …) — it falls through to the random-sampling stage.
  There, per trial it draws `sample_random_value()`s, `substitute_symbols` →
  `sub`, then `decide_numeric(sub)` numericalizes `sub` at each precision rung.
  Some intermediate `Expr` tree built while numericalizing a not-fully-numeric
  sample is **not freed** — the leaked roots are small `Plus/Times` trees that
  still contain symbols (i.e. partially-substituted / partially-numericalized
  subexpressions), pointing at `numericalize`/`evaluate_rung` rather than at
  `sub`/`vals` (both of which `decide_schwartz_zippel` does free correctly at
  lines 686/689).
- **Trigger / nondeterminism:** the leak depends on which random values are
  drawn (`sample_random_value` sometimes returns a Real, sometimes a Complex,
  exercising different numericalize branches), so it is **intermittent** — a
  single call may leak 0; many calls reliably accumulate. Independent of
  `Integrate`: reproduced directly with bare `PossibleZeroQ` calls:
  ```bash
  { for i in $(seq 1 40); do
      echo "PossibleZeroQ[(x+$i)^(1/3) + Sqrt[x+$i] - x^(1/2)]"; done; echo 'Quit[]'; } \
    | leaks --atExit -- ./Mathilda 2>&1 | grep 'leaks for'
  # => ~1300 leaks / ~70 KB (vs 0 for the `1+1` baseline)
  ```
- **Stack (representative, inlined):**
  ```
  malloc → expr_new_function (src/expr.c:81)
  zero_test_decide (src/zero_test.c:720)        # inlined decide_schwartz_zippel→decide_numeric→numericalize
  builtin_possible_zero_q (src/zero_test.c:734)
  evaluate_step (src/eval.c:813)
  evaluate (src/eval.c:1079)
  ...
  ```
- **Likely fix:** audit `numericalize` (`src/numeric.c`) for an intermediate
  `Expr` it allocates but doesn't free when the input is only partially numeric
  (head/argument rebuild path). The SZ driver itself (`decide_schwartz_zippel`)
  already balances `sub` and `vals`; the unfreed tree is created one level
  deeper, inside `numericalize_at(sub, bits)`. A scoped reproduction is to call
  `numericalize` on `Plus[Real, Power[x, Rational[1,3]]]` (a Plus mixing a Real
  and a symbolic radical) and leak-check that in isolation.
- **Severity:** low–moderate. Small per leak (~48–320 B per root tree) but
  **unbounded in aggregate** — `PossibleZeroQ` is on hot verification paths
  (every `Integrate` differentiate-back gate, `Simplify`, etc.), so a session
  that does much symbolic-radical integration accumulates KBs.
- **Discovered:** 2026-06-06, surfaced by `Integrate`LinearRadicals` — its
  differentiate-back verification (identical to `DerivativeDivides`) feeds
  radical `D[result,x] - f` residues to `PossibleZeroQ`; the DerivativeDivides
  test integrands happen to settle before Stage 3, so they never exposed it.
