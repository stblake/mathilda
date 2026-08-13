# Head

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Head[expr]`**

gives the head of expr.

**`Head[expr, h]`**

wraps the result with h, i.e. returns h\[Head\[expr\]\].

<details>
<summary>Notes</summary>

For atoms, Head returns Integer, Real, BigInt, Rational, Complex, Symbol, or String; for a compound expression f\[...\], Head returns f.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Head[f[x]]
Out[1]= f

In[2]:= Head[3/4]
Out[2]= Rational

In[3]:= Head[a + b, f]
Out[3]= f[Plus]
```

## Implementation notes

`builtin_head` (in `src/part.c`) returns the head of its argument via `expr_head` — `f` for `f[...]`, and the type symbol (`Integer`, `Symbol`, `List`, etc.) for atoms. The 2-arg form `Head[expr, h]` wraps the result as `h[Head[expr]]`, leaving the outer application for the evaluator to reduce.

- For functions, returns the symbol or expression acting as the head.
- For atoms, returns the symbolic type name: `Integer`, `Real`, `Rational`, `Complex`, `Symbol`, or `String`.

**Attributes:** none registered.

## References

**See also:** [Rational](../../arithmetic/Rational/), [Complex](../../arithmetic/Complex/), [Symbol](../../expression-information/Symbol/)

- Source: [`src/part.c`](https://github.com/stblake/mathilda/blob/main/src/part.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_blas.c`](https://github.com/stblake/mathilda/blob/main/tests/test_blas.c)
- Tests: [`tests/test_cherry_dilog.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_dilog.c)
