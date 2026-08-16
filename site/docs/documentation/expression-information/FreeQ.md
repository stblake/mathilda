# FreeQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FreeQ[expr, form]`**

yields True if no subexpression of expr matches form, False otherwise.

**`FreeQ[expr, form, levelspec]`**

restricts the search to parts of expr at the levels specified by levelspec.

**`FreeQ[form]`**

is the operator form: FreeQ\[form\]\[expr\] == FreeQ\[expr, form\].

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= FreeQ[{1, 2, 4, 1, 0}, 0]
Out[1]= False

In[2]:= FreeQ[{a, b, b, a, a, a}, _Integer]
Out[2]= True

In[3]:= f[c_ x_, x_] := c f[x, x] /; FreeQ[c, x]

In[4]:= {f[3 x, x], f[a x, x], f[(1 + x) x, x]}
Out[4]= {3 f[x, x], a f[x, x], f[x (1 + x), x]}
```

### Applications (4)

```mathematica
In[5]:= FreeQ[x^2 + y^2, z]
Out[5]= True

In[6]:= FreeQ[x^2 + y^2, y]
Out[6]= False

In[7]:= FreeQ[D[Sin[x] Exp[x], x], Cos]
Out[7]= False

In[8]:= Select[Range[20], FreeQ[FactorInteger[#], {2, _}] &]
Out[8]= {1, 3, 5, 7, 9, 11, 13, 15, 17, 19}
```

## Implementation notes

`builtin_freeq` (`src/funcprog.c`) returns `True` iff no subexpression of the first argument matches the pattern `form`. `freeq_at_level` walks the tree depth-first, calling the pattern matcher `match` at each level permitted by the level spec (default `{0, Infinity}`), and short-circuits to `False` on the first match. `Rational`/`Complex` are treated as atomic; `Heads -> False` excludes function heads from the search.

- `Protected`.
- By default, explores levels `{0, Infinity}` and option `Heads -> True` is enabled.
- `form` can be a structural pattern.

**Attributes:** `Protected`.

## References

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_cherry_dilog.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_dilog.c)
- Tests: [`tests/test_cherry_ei.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_ei.c)
- Tests: [`tests/test_cherry_li.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_li.c)
- Tests: [`tests/test_condition_downvalue.c`](https://github.com/stblake/mathilda/blob/main/tests/test_condition_downvalue.c)

## Notes & additional examples

### Notes

`FreeQ[expr, form]` tests whether *no* subexpression of `expr` matches the
pattern `form`, searching at every level. `form` may be a literal symbol or
a full pattern: the third example confirms that differentiating
`Sin[x] Exp[x]` introduces a `Cos` term. The last example is a structural
sieve — keeping only integers whose `FactorInteger` output is free of any
`{2, _}` pair recovers the odd numbers, with `FreeQ` matching against the
`{prime, exponent}` factor structure rather than a numeric value.
