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

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= r = 0; Do[If[EvenQ[i], Continue[]]; r += i, {i, 10}]; r
Out[1]= 25

In[2]:= r = 0; For[i = 1, i <= 10, i++, If[EvenQ[i], Continue[]]; r += i]; r
Out[2]= 25
```

## Implementation notes

- Has attribute `Protected`.
- Takes effect as soon as it is evaluated. In `Do` it advances the iterator and
  re-tests; in `For` it evaluates the increment step then re-tests; in `While` it
  re-evaluates the test.
- Outside any loop, `Continue[]` emits the message `Continue::nofwd` and returns
  `Hold[Continue[]]`.

**Attributes:** `Protected`.

## References

**See also:** [Do](../../control-flow/Do/), [For](../../control-flow/For/), [While](../../control-flow/While/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_iter.c`](https://github.com/stblake/mathilda/blob/main/tests/test_iter.c)
- Tests: [`tests/test_return.c`](https://github.com/stblake/mathilda/blob/main/tests/test_return.c)
