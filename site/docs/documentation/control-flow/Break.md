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

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Do[Print[i]; If[i > 2, Break[]], {i, 10}] 1 2 3
Out[1]= 6 Null

In[2]:= For[i = 1, i <= 10, i++, If[i > 2, Break[]]]; i
Out[2]= 3
```

## Implementation notes

- Has attribute `Protected`.
- Takes effect as soon as it is evaluated (e.g. inside an `If` within the body),
  escaping only the *innermost* enclosing loop.
- Outside any loop, `Break[]` emits the message `Break::nofwd` and returns
  `Hold[Break[]]` (inert, so feeding it back does not re-trigger).

**Attributes:** `Protected`.

## See also

[Do](../../control-flow/Do/), [For](../../control-flow/For/), [While](../../control-flow/While/), [If](../../control-flow/If/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_iter.c`](https://github.com/stblake/mathilda/blob/main/tests/test_iter.c)
