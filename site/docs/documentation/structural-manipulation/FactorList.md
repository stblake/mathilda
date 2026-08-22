# FactorList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FactorList[poly] gives a list of the irreducible factors of poly together`**

<details>
<summary>Notes</summary>

with their exponents, as {factor, exponent} pairs.  A thin wrapper over Factor: options (GaussianIntegers -\> True, Extension -\> {a1, ...}) are forwarded verbatim.  The first element is always the overall numerical factor {c, 1} (it is {1, 1} when there is none); denominator factors of a rational function appear with negative exponents.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= FactorList[x^2 - 1]
Out[1]= {{1, 1}, {-1 + x, 1}, {1 + x, 1}}

In[2]:= FactorList[2 x^3 + 2 x^2 - 2 x - 2]
Out[2]= {{2, 1}, {-1 + x, 1}, {1 + x, 2}}

In[3]:= FactorList[(x^3 + 2 x^2)/(x^2 - 4 y^2) - (x + 2)/(x^2 - 4 y^2)]
Out[3]= {{1, 1}, {-1 + x, 1}, {1 + x, 1}, {2 + x, 1}, {x - 2 y, -1}, {x + 2 y, -1}}
```

### Options (1)

```mathematica
In[4]:= FactorList[x^2 + 1, Extension -> I]
Out[4]= {{1, 1}, {-I + x, 1}, {I + x, 1}}
```

## Algorithm

facpoly_list.c -- FactorList[poly] / FactorList[poly, opts].

A thin wrapper over Factor: factor `poly` via `Factor[poly, opts...]` (all options are forwarded verbatim -- GaussianIntegers, Extension, ...), then split the resulting product into {factor, exponent} pairs.

```text
  FactorList[x^2 - 1]              -> {{1, 1}, {-1 + x, 1}, {1 + x, 1}}
  FactorList[2 x^3 + 2 x^2 - ...]  -> {{2, 1}, {-1 + x, 1}, {1 + x, 2}}
  FactorList[x^4 - 2, Extension -> Sqrt[2]]
                                   -> {{1, 1}, {Sqrt[2] + x^2, 1}, {-Sqrt[2] + x^2, 1}}
```

The first element is always the overall numerical factor {c, 1} (which is

```text
{1, 1} when there is no numerical factor).  Denominator factors of a
```

rational function carry negative exponents.

Parsing rules on the Factor output R (a Times, a bare factor, or a number):

```text
  - a number literal (Integer / Rational / Real / Complex / ...) multiplies
    into the overall numerical factor `c`;
  - Power[base, e] with an *integer* e is the pair {base, e} (this is a
    factor raised to a multiplicity, positive or negative);
  - anything else -- including Power[base, 1/2] = Sqrt[base], which is an
    irreducible factor in its own right -- is the pair {factor, 1}.
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Factor x^120 - 1 | 7.01 s | 0.044 s | 3.85 s |
| Factor dense degree 60 | 2.2 s | 0.599 s | 26.5 s |
| Factor product of 16 quadratics | 1.29 s | 0.722 s | 27.6 s |
| Factor sparse degree 60 | 0.532 s | 0.516 s | 9.48 s |
| Factor product of 8 quadratics | 0.46 s | 0.191 s | 3.31 s |
| Factor 5-var, 4 factors | 0.377 s | 0.415 s | 20 s |

## Implementation notes

- `Listable`, `Protected`.
- A thin wrapper over `Factor`: it factors via `Factor[poly, opts...]` (options are forwarded verbatim) and splits the product into `{factor, exponent}` pairs.
- The first element is always the overall numerical factor `{c, 1}` — it is `{1, 1}` when there is no numerical factor.
- Denominator factors of a rational function appear with negative exponents.
- `Times @@ Power @@@ FactorList[poly]` reconstructs `poly` (up to `Factor`'s normal form).
- Because it delegates to `Factor`, the factorisation depth and normal form (factor ordering, `Extension`/`GaussianIntegers` handling) are exactly `Factor`'s.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [Factor](../../algebra/Factor/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_factorlist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_factorlist.c)
