# Exponent

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Exponent[expr, form] gives the maximum power with which form appears in the`**

**`Exponent[expr, form, h] applies h to the set of exponents (default Max).`**

**`Exponent[0, x] is -Infinity.  Exponent[expr, {f1, f2, ...}] gives the list of`**

<details>
<summary>Notes</summary>

expanded form of expr. form may be a symbol, a kernel, or a product of terms; expr need not be expanded.  Exponent is purely syntactic (no zero-coefficient recognition). exponents for each fi.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Exponent[1 + x^2 + a x^3, x]
Out[1]= 3

In[2]:= Exponent[0, x]
Out[2]= -Infinity

In[3]:= Exponent[x^(n+1) + 2 Sqrt[x] + 1, x]
Out[3]= Max[1/2, 1 + n]

In[4]:= Exponent[(x^2+1)^3 - 1, x, Min]
Out[4]= 2

In[5]:= Exponent[1 + x^2 + a x^3, x, List]
Out[5]= {0, 2, 3}
```

## Algorithm

exponent.c -- Exponent[expr, form] / Exponent[expr, form, h].

Gives the maximum power (default h = Max) with which `form` appears in the

```text
expanded form of `expr`.  Purely syntactic: `expr` is expanded first, but no
```

zero-coefficient recognition is attempted (Exponent[zero x^2 + x + 1, x] = 2 even when `zero` is numerically 0 but not in normal form).

Algorithm ---------

```text
  1. expanded = Expand[expr].
  2. Split into additive terms (the args of a top-level Plus, else the whole
     expression as a single term).  The genuine zero polynomial has NO terms,
     so its exponent set is empty and h[] fires -- Max[] = -Infinity.
  3. For each term (a monomial), read off the exponent of `form`:
       - form decomposes into (base, fe) pairs (a symbol/kernel -> (form,1);
         a Power[b,e] -> (b,e); a product -> one pair per factor);
       - the exponent of a single base in a monomial is the power to which it
         is raised (0 if absent, symbolic/rational allowed);
       - for a product form the term's exponent is Min over the form's bases
         of (base-exponent / fe) -- the largest k with form^k dividing it.
  4. Collect the exponents into a sorted, de-duplicated set and return
     h @@ set (h defaults to Max).  Symbolic exponents flow through Max/Min
     unevaluated (Max[1/2, 1 + n]).

Listable + Protected.  Listable makes Exponent[expr, {f1, f2, ...}] thread
```

into the per-form list of exponents for free.

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Discriminant of deg 20 | 2.51 s | 0.068 s | 0.182 s |
| Expand (1+x)^400 | 0.434 s | 0.107 s | 0.003 s |
| Cancel deg-60 over deg-58 | 0.337 s | 0.569 s | 7.37 s |
| PolynomialGCD, coprime deg 40 | 0.252 s | 0.087 s | 0.334 s |
| PolynomialGCD, shared deg-20 factor | 0.079 s | 0.063 s | 0.764 s |
| PolynomialQuotient deg 60 / deg 20 | 0.063 s | 0.209 s | 0.945 s |

## Implementation notes

- `Listable`, `Protected`.
- The default aggregator is `h = Max`. `Exponent[expr, form, Min]` gives the lowest power; `Exponent[expr, form, List]` gives the sorted, de-duplicated set of exponents.
- `form` may be a symbol, a kernel (e.g. `Sin[x]`), or a product of terms.
- Works whether or not `expr` is explicitly given in expanded form (it expands internally).
- Purely syntactic: it does not attempt to recognise zero coefficients.
- Exponents may be rational numbers or symbolic expressions.
- `Exponent[0, x]` is `-Infinity` (empty exponent set, `h = Max`).
- The `Listable` attribute makes `Exponent[expr, {form1, form2, ...}]` give the list of exponents for each `formi`.

**Attributes:** `Listable`, `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_characteristicpolynomial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_characteristicpolynomial.c)
- Tests: [`tests/test_expand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expand.c)
- Tests: [`tests/test_exponent.c`](https://github.com/stblake/mathilda/blob/main/tests/test_exponent.c)
- Tests: [`tests/test_subresultantpolynomials.c`](https://github.com/stblake/mathilda/blob/main/tests/test_subresultantpolynomials.c)
