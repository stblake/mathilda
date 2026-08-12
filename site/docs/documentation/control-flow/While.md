# While

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`While[test, body] evaluates test, then body, repeatedly, until test first fails to give True.`**

**`While[test] does the loop with a Null body, which is useful when test has side-effects.`**

**`Break[] inside body exits the loop.`**

**`Continue[] inside body skips the rest of body and re-evaluates test.`**

**`Return[v] inside body causes While to yield v; otherwise While returns Null.`**

<details>
<summary>Notes</summary>

While has attribute HoldAll.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= n = 1; While[n < 4, n = n + 1]; n
Out[1]= 4

In[2]:= {a, b} = {27, 6}; While[b != 0, {t1, t2} = {b, Mod[a, b]}; a = t1; b = t2]; a
Out[2]= 3

In[3]:= n = 1; While[True, If[n > 10, Break[]]; n = n + 1]; n
Out[3]= 11
```

## Implementation notes

**Algorithm.** `builtin_while` implements `While[test]` or `While[test, body]` and is `ATTR_HOLDALL`, so both arguments are re-evaluated each iteration. Each pass calls `evaluate(test)`; the result is first run through `iter_flow_classify` (boundary head `SYM_While`) so flow-control raised inside the test itself is honored (Break/Return exit, Throw/Abort/Quit propagate, Continue restarts the loop). The loop then exits unless the test is exactly the symbol `True`. The body (if present) is evaluated and classified the same way: `Break` exits (yielding `Null`), `Continue` skips to the next test, `Return[v]` exits yielding `v`. Truth testing is pointer equality on `SYM_True`. The loop returns the Return payload if one was issued, else `Null`.

- Has attribute `HoldAll`; both `test` and `body` are re-evaluated each iteration.
- `Break[]` inside `body` exits the loop, yielding `Null`.
- `Continue[]` inside `body` skips the rest of `body` and returns to re-evaluating `test`.
- `Return[v]` inside `body` causes `While` to yield `v`.
- `Throw`, `Abort`, and `Quit` propagate unchanged.
- If the very first evaluation of `test` is not `True`, `body` is never evaluated.
- Returns `Null` unless an explicit `Return` is issued.

**Attributes:** `HoldAll`, `Protected`.

## See also

[HoldAll](../../expression-information/HoldAll/), [Throw](../../control-flow/Throw/), [Return](../../control-flow/Return/)

## References

- Source: [`src/iter.c`](https://github.com/stblake/mathilda/blob/main/src/iter.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_iter.c`](https://github.com/stblake/mathilda/blob/main/tests/test_iter.c)
- Tests: [`tests/test_mateigen_direct.c`](https://github.com/stblake/mathilda/blob/main/tests/test_mateigen_direct.c)
