# Identity

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Identity[expr] gives expr unchanged (the identity function).`**

## Examples (6)

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

### Applications (3)

Returns its argument unchanged

```mathematica
In[4]:= Identity[x]
Out[4]= x
```

Its argument is evaluated normally first

```mathematica
In[5]:= Identity[1 + 1]
Out[5]= 2
```

The identity function is handy as a default callback

```mathematica
In[6]:= Map[Identity, {1, 2, 3}]
Out[6]= {1, 2, 3}
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

## Notes & additional examples

### Notes

`Identity[expr]` returns `expr` unchanged. It takes exactly one argument; any other arity
is left unevaluated. Having no held attributes, it evaluates its argument through the
normal pipeline before returning it, so `Identity[1 + 1]` is `2`.

It is most useful as a neutral function argument — a do-nothing callback where some
transformation is expected.
