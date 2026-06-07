# Symbol

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Symbol["name"]
    refers to the symbol with the specified name, creating it in
    $Context if none yet exists.

All symbols, whether explicitly entered using Symbol or not, have head
Symbol; x_Symbol matches any symbol. The name string may contain
letters, letter-like forms, or digits but must not start with a digit.
A backtick (`) separates context prefixes; a leading backtick makes the
name relative to the current context $Context.

Attributes: Protected.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Symbol["x"]
Out[1]= x

In[2]:= Head[%]
Out[2]= Symbol

In[3]:= {f[x], f["x"], f[2]} /. f[s_Symbol] :> g[s]
Out[3]= {g[x], f["x"], f[2]}

In[4]:= Symbol["a`x"]
Out[4]= a`x

In[5]:= Symbol["1x"]
Out[5]= Symbol["1x"]
```

## Implementation notes

- `Protected`.
- Every expression's `Head` matches `Symbol` for symbols; `x_Symbol` patterns therefore match any symbol.
- The string must satisfy the standard symbol-name syntax: each segment (separated by backticks) starts with a letter or `$`, followed by letters, digits, or `$`.
- A leading backtick (`Symbol["\`x"]`) makes the name relative to the current `$Context`. An embedded backtick (`Symbol["a\`x"]`) is treated as an absolutely-qualified name. A bare name is resolved through the standard `$Context` / `$ContextPath` rules.
- Invalid names emit `Symbol::symname` to `stderr` and leave the call unevaluated; non-string arguments also leave the call unevaluated.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
