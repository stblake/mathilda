# Break

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Break[] exits the nearest enclosing Do, For, or While loop.`**

**`Break[] takes effect as soon as it is evaluated.`**

<details>
<summary>Notes</summary>

After Break\[\], the enclosing loop returns Null. Break has attribute Protected.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Do[Print[i]; If[i > 2, Break[]], {i, 10}] 1 2 3
Out[1]= 6 Null

In[2]:= For[i = 1, i <= 10, i++, If[i > 2, Break[]]]; i
Out[2]= 3
```

### Applications (4)

The loop stops the moment the guard fires

```mathematica
In[3]:= s = 0; Do[If[i > 3, Break[]]; s = s + i, {i, 10}]; s
Out[3]= 6
```

Break escapes the loop, leaving i at its exit value

```mathematica
In[4]:= For[i = 1, i <= 10, i++, If[i == 4, Break[]]]; i
Out[4]= 4
```

The standard way out of a While[True] loop

```mathematica
In[5]:= n = 1; While[True, If[n > 5, Break[]]; n++]; n
Out[5]= 6
```

Outside any loop it is inert: reported and wrapped in Hold

```mathematica
In[6]:= Break[]
Out[6]= Hold[Break[]]
```

## Implementation notes

**Algorithm.** `Break[]` is a zero-argument flow-control marker, `Protected` with
no `Hold` attributes. `builtin_break` only validates arity — a non-zero argument
count emits the standard `Break::argx` message through `builtin_arg_error` and
leaves the call unevaluated — and otherwise returns `NULL`, so the raw `Break[]`
node stands unevaluated as the marker itself. It belongs to the **head-detected**
family of non-local control (Mechanism B: `Return`/`Break`/`Continue`/`Abort`/
`Quit`), distinct from the `Throw`/`Goto` sentinel. The difference is where it is
seen: `Do`/`For`/`While` (`src/iter.c`) inspect it by head through
`iter_flow_classify` (keyed on the interned `SYM_Break`) at the loop boundary,
but it is **not** short-circuited in `evaluate_step`'s argument-evaluation loop,
so `Print[Break[]]` does not escape — the semantics are lexical.

**Data structures.** None; `Break[]` carries no payload and is just an unevaluated
`EXPR_FUNCTION` node that the loops match by pointer-identity of its head symbol.

**Complexity / limits.** A `Break[]` takes effect as soon as it is evaluated and
escapes only the *innermost* enclosing `Do`/`For`/`While`, which then yields
`Null`. `Table` deliberately does not honour it. A `Break[]` that reaches top
level with no loop to consume it is reported by
`eval_report_uncaught_break_continue` (`src/eval.c`) with `Break::nofwd` and
rewritten to the inert `Hold[Break[]]`, so feeding it back does not re-trigger.

- Has attribute `Protected`.
- Takes effect as soon as it is evaluated (e.g. inside an `If` within the body),
  escaping only the *innermost* enclosing loop.
- Outside any loop, `Break[]` emits the message `Break::nofwd` and returns
  `Hold[Break[]]` (inert, so feeding it back does not re-trigger).

**Attributes:** `Protected`.

## References

**See also:** [Do](../../control-flow/Do/), [For](../../control-flow/For/), [While](../../control-flow/While/), [If](../../control-flow/If/)

- Source: [`src/iter.c`](https://github.com/stblake/mathilda/blob/main/src/iter.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_iter.c`](https://github.com/stblake/mathilda/blob/main/tests/test_iter.c)

## Notes & additional examples

### Notes

`Break[]` takes effect as soon as it is evaluated — even from inside an `If`
within the body — and escapes only the *innermost* enclosing `Do`, `For` or
`While`, which then yields `Null`.

It is recognised by its head at the loop boundary, not short-circuited the way a
`Throw` is, so `Print[Break[]]` does **not** escape the loop. `Table` deliberately
does not honour `Break`. A `Break[]` that reaches top level with no loop to
consume it prints `Break::nofwd` and comes back as the inert `Hold[Break[]]`, so
feeding that result back in does not re-trigger anything.
