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
