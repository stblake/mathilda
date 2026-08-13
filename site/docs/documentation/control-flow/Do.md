# Do

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Do[expr, n] evaluates expr n times.`**

**`Do[expr, {i, imax}] evaluates expr with i successively taking on values 1 through imax.`**

**`Do[expr, {i, imin, imax}] starts with i = imin.`**

**`Do[expr, {i, imin, imax, di}] uses steps di.`**

**`Do[expr, {i, {i1, i2, ...}}] uses the successive values i1, i2, ....`**

**`Do[expr, {n}] evaluates expr n times with no iteration variable.`**

**`Do[expr, iter1, iter2, ...] iterates over multiple variables, with the rightmost varying fastest.`**

**`Break[] inside expr exits the innermost Do loop.`**

**`Continue[] inside expr skips the rest of expr and proceeds to the next iteration.`**

**`Return[v] inside expr causes the enclosing function to yield v; Do itself returns Null.`**

<details>
<summary>Notes</summary>

Do has attribute HoldAll: expr and the iterator specifications are held unevaluated until each iteration.

</details>

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

```mathematica
In[1]:= (s = 0; Do[s = s + i, {i, 5}]; s)
Out[1]= 15

In[2]:= (p = 1; Do[p = p*2, {3}]; p)
Out[2]= 8
```

## Implementation notes

**Algorithm.** `builtin_do` is `ATTR_HOLDALL`, so the body and iterator spec are re-evaluated on every pass rather than once at call time. A multi-spec `Do[expr, s1,...,sk]` is rewritten as nested two-spec `Do[Do[expr, sk], s1,...,s_{k-1}]` and handed back to the evaluator, reducing every case to a single iterator. The spec is parsed by the shared `iter_spec_parse` into an `IterSpec` (kinds `COUNT` for `n`/`{n}`, `RANGE` for `{i,imin,imax,di}` with defaults, `LIST` for `{i,{...}}`); numeric bounds are resolved to doubles by `iter_spec_resolve_numeric`, with `Infinity` allowed for unbounded loops.

**Variable localization.** Before looping, `iter_spec_shadow` saves and clears the iterator symbol's `own_values`; each iteration binds the current value via `symtab_add_own_value`, and `iter_spec_restore` frees the per-iteration binding chain and restores the original OwnValue afterward — a manual mimic of standard iterator scoping. For exact ranges the bound value `curr_e` is advanced with `Plus[curr_e, di_e]` evaluated each step (keeping integer/rational exactness), while a parallel `double val` drives the termination comparison (with a `1e-14` tolerance).

**Control flow.** After each `evaluate(body)`, `iter_flow_classify` (boundary head `SYM_Do`) maps the result to break / continue / return-value / propagate (Throw/Abort/Quit/foreign Return). `Continue` in a range loop still advances the counter before re-testing. Returns the Return payload if any, else `Null`.

- `HoldAll`, evaluating its body only after arguments are substituted.
- Employs exact dynamic iteration identical to `Table` but discards the evaluated results, returning `Null`.
- Supports explicit break states (`Return`, `Break`, `Continue`, `Throw`, `Abort`, `Quit`).
- Can execute an infinite loop using `Do[expr, Infinity]`.
- A body that is machine-numeric throughout takes an automatic fast path
  (`src/numloop.c`) that runs it as a double-stack program with no `Expr`
  allocation. The fast path is built without evaluating anything: loop-invariant
  terms are folded only when they are syntactically numeric (literals, `Pi`,
  `Sqrt[2]`, `p/q`, arithmetic and elementary heads). A body containing anything
  else — notably a call to a user-defined function — declines the fast path and
  runs interpreted, so `expr` is evaluated exactly as many times as the iterator
  specifies and its side effects fire exactly that often. The same rule governs
  `For`, `While`, `Nest`, `NestWhile`, `FixedPoint`, `Fold` and `Map`.
- An **exact-integer** counter body (`Do[s = s + i, {i, 1, n}]`, `For[...]`, no
  inexact leaf) takes the same fast path in overflow-checked `int64` and returns
  an exact `Integer` — the double stack is used only when the body is inexact. On
  an `int64` overflow, or a step that is not integer-closed (a `Rational` from a
  division, a transcendental), the whole run bails to the interpreter, so the
  answer is always the interpreter's: `Do[p = p*i, {i, 1, 25}]` is the full `25!`
  bignum, `Do[s = s + i/2, {i, 1, 10}]` is `55/2`, never a wrapped machine
  integer.
- An **integer range** terminates on an exact comparison of the running value
  against the bound, so a span anywhere in the `int64` range is correct — including
  at the very top, where consecutive values differ by less than the `double`
  spacing: `Do[…, {i, 9223372036854775805, 9223372036854775807}]` runs exactly 3
  times. This holds on all three paths — interpreter, the `numloop.c` fast path
  (whose counter increment is overflow-checked so it stops at the edge rather than
  wrapping), and `Compile[]` (whose loop-control arithmetic stays overflow-checked
  even under `RuntimeOptions -> "Speed"`, bailing to the interpreter at the
  boundary; only *value* arithmetic wraps in that mode). `Sum`, `Product` and
  `Table` share the same termination.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [Table](../../lists-and-iteration/Table/), [Return](../../control-flow/Return/), [Break](../../control-flow/Break/), [Continue](../../control-flow/Continue/), [Throw](../../control-flow/Throw/), [Pi](../../mathematical-constants/Pi/), [For](../../control-flow/For/)

- Source: [`src/iter.c`](https://github.com/stblake/mathilda/blob/main/src/iter.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_catch_throw.c`](https://github.com/stblake/mathilda/blob/main/tests/test_catch_throw.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)

## Notes & additional examples

### Notes

`Do[expr, {i, imax}]` runs `expr` with `i` taking values `1` through `imax`; `Do[expr, n]` simply repeats `expr` `n` times. `Do` returns `Null`, so it is used for side effects — read the accumulated value out of a variable afterwards.
