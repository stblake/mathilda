# Continue

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Continue[] proceeds to the next iteration of the nearest enclosing Do, For, or While loop.`**

**`Continue[] skips the remainder of the current loop body.`**

**`Continue[] takes effect as soon as it is evaluated.`**

<details>
<summary>Notes</summary>

Continue has attribute Protected.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= r = 0; Do[If[EvenQ[i], Continue[]]; r += i, {i, 10}]; r
Out[1]= 25

In[2]:= r = 0; For[i = 1, i <= 10, i++, If[EvenQ[i], Continue[]]; r += i]; r
Out[2]= 25
```

### Applications (3)

Skip the even i and sum the odd ones

```mathematica
In[3]:= r = 0; Do[If[EvenQ[i], Continue[]]; r += i, {i, 10}]; r
Out[3]= 25
```

In For, Continue still runs the increment step

```mathematica
In[4]:= r = 0; For[i = 1, i <= 6, i++, If[i == 3, Continue[]]; r += i]; r
Out[4]= 18
```

Inert outside any loop

```mathematica
In[5]:= Continue[]
Out[5]= Hold[Continue[]]
```

## Implementation notes

**Algorithm.** `Continue[]` is a zero-argument flow-control marker, `Protected`
with no `Hold` attributes. `builtin_continue` validates arity only — a non-zero
argument count emits `Continue::argx` through `builtin_arg_error` — and otherwise
returns `NULL`, so the raw `Continue[]` node stands as the marker. Like `Break`,
it is a **head-detected** marker (Mechanism B), recognised by the loop builtins
through `iter_flow_classify` (interned `SYM_Continue`) and not short-circuited in
the argument-evaluation loop, so it only takes effect at a loop boundary.

**Data structures.** None; a bare unevaluated node matched by head.

**Complexity / limits.** `Continue[]` skips the remainder of the current body of
the innermost `Do`/`For`/`While` and advances to the next iteration: in `Do` it
advances the iterator and re-tests (the arithmetic-progression form must still
step its running value on `ITER_FLOW_CONTINUE`, else the loop would spin), in
`For` it evaluates the increment step and re-tests, and in `While` it
re-evaluates the test. A `Continue[]` that reaches top level is reported with
`Continue::nofwd` by `eval_report_uncaught_break_continue` and rewritten to the
inert `Hold[Continue[]]`.

- Has attribute `Protected`.
- Takes effect as soon as it is evaluated. In `Do` it advances the iterator and
  re-tests; in `For` it evaluates the increment step then re-tests; in `While` it
  re-evaluates the test.
- Outside any loop, `Continue[]` emits the message `Continue::nofwd` and returns
  `Hold[Continue[]]`.

**Attributes:** `Protected`.

## References

**See also:** [Do](../../control-flow/Do/), [For](../../control-flow/For/), [While](../../control-flow/While/)

- Source: [`src/iter.c`](https://github.com/stblake/mathilda/blob/main/src/iter.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_iter.c`](https://github.com/stblake/mathilda/blob/main/tests/test_iter.c)
- Tests: [`tests/test_return.c`](https://github.com/stblake/mathilda/blob/main/tests/test_return.c)

## Notes & additional examples

### Notes

`Continue[]` abandons the rest of the current loop body and moves to the next
iteration of the innermost `Do`, `For` or `While`. What "next iteration" means
differs by loop: `Do` advances its iterator and re-tests, `For` evaluates its
increment step and re-tests, and `While` re-evaluates its test.

Like `Break`, it is a head-detected marker and takes effect only at a loop
boundary. Outside any loop it prints `Continue::nofwd` and returns the inert
`Hold[Continue[]]`.
