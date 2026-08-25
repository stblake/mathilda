# Symbol

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Symbol["name"]`**

refers to the symbol with the specified name, creating it in $Context if none yet exists.

<details>
<summary>Notes</summary>

All symbols, whether explicitly entered using Symbol or not, have head Symbol; x\_Symbol matches any symbol. The name string may contain letters, letter-like forms, or digits but must not start with a digit. A backtick (\`) separates context prefixes; a leading backtick makes the name relative to the current context $Context. Attributes: Protected.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Symbol["x"]
Out[1]= x

In[2]:= Head[%]
Out[2]= Out

In[3]:= {f[x], f["x"], f[2]} /. f[s_Symbol] :> g[s]
Out[3]= {g[x], f["x"], f[2]}

In[4]:= Symbol["a`x"]
Out[4]= a`x
```

## Implementation notes

`builtin_symbol` (`src/core.c`) converts a string to a symbol. It validates the name with `symbol_name_is_valid` (each backtick-delimited context segment must start with a letter or `$` and continue with alphanumerics/`$`), runs it through `context_resolve_name` to apply the current context, and returns an `EXPR_SYMBOL`. Invalid names emit a `Symbol::symname` message and return `NULL`.

- `Protected`.
- Every expression's `Head` matches `Symbol` for symbols; `x_Symbol` patterns therefore match any symbol.
- The string must satisfy the standard symbol-name syntax: each segment (separated by backticks) starts with a letter or `$`, followed by letters, digits, or `$`.
- A leading backtick (`Symbol["\`x"]`) makes the name relative to the current `$Context`. An embedded backtick (`Symbol["a\`x"]`) is treated as an absolutely-qualified name. A bare name is resolved through the standard `$Context` / `$ContextPath` rules.
- Invalid names emit `Symbol::symname` to `stderr` and leave the call unevaluated; non-string arguments also leave the call unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [Head](../../structural-manipulation/Head/), [$Context](../../scoping-constructs/$Context/), [$ContextPath](../../scoping-constructs/$ContextPath/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_findmin_dogleg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findmin_dogleg.c)
- Tests: [`tests/test_findmin_methods.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findmin_methods.c)
- Tests: [`tests/test_findmin_neldermead.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findmin_neldermead.c)
- Tests: [`tests/test_findmin_newtoncg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findmin_newtoncg.c)
