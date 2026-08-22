# Identity

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Identity[expr] gives expr unchanged (the identity function).`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Identity[x]
Out[1]= x

In[2]:= Identity[1 + 2]
Out[2]= 3

In[3]:= Map[Identity, {a, b, c}]
Out[3]= {a, b, c}
```

## Implementation notes

`builtin_identity` (`src/core.c`) is the one-argument identity: it returns a copy of its single argument unchanged, or `NULL` for any other arity.

**Attributes:** `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_eval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eval.c)
- Tests: [`tests/test_repl_hooks.c`](https://github.com/stblake/mathilda/blob/main/tests/test_repl_hooks.c)
- Tests: [`tests/test_sequence.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sequence.c)
