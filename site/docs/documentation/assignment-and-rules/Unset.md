# Unset

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Unset[lhs] or lhs =.`**

removes any rule whose left-hand side is lhs, up to renaming of pattern variables. A bare symbol clears its value; a function form clears the matching definition on the head symbol.

<details>
<summary>Notes</summary>

Unset has attribute HoldFirst; Protected symbols are not affected.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= x = 5; x =.; x
Out[1]= x

In[2]:= f[x_] := x^2; f[x_] =.; f[3]
Out[2]= f[3]

In[3]:= fact[1] = 1; fact[n_] := n fact[n - 1]; fact[1] =.; fact[1]
Out[3]= 0
```

## Implementation notes

- `=.` is a low-precedence postfix operator (precedence 40, like `Set`), so it
  captures the whole preceding expression: `a b =.` parses as `Unset[a b]`. The
  guard against a trailing digit keeps `k =.5` parsing as `Set[k, 0.5]`.
- `Unset` has attributes `{HoldFirst, Protected}`; it holds `lhs`, so the symbol
  (not its value) is operated on. `Protected`/`Locked` symbols are not affected.
- Always returns `Null`, whether or not a matching rule was found.

**Attributes:** `HoldFirst`, `Protected`.

## References

**See also:** [Set](../../assignment-and-rules/Set/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_unset.c`](https://github.com/stblake/mathilda/blob/main/tests/test_unset.c)
