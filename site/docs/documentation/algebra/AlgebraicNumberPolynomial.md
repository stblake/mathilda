# AlgebraicNumberPolynomial

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AlgebraicNumberPolynomial[a, x]`**

gives the polynomial in x corresponding to the AlgebraicNumber object a.

<details>
<summary>Notes</summary>

For a = AlgebraicNumber\[theta, {c0, c1, ..., cn}\], the result is the polynomial c0 + c1 x + ... + cn x^n, from which a is recovered by replacing x with theta.  An integer or rational a is the constant polynomial and is returned unchanged; any other argument stays unevaluated.  Threads over lists.  Attributes: Listable, Protected.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2] + Sqrt[3], {1, 2, 3, 4}], x]
Out[1]= 1 + 2 x + 3 x^2 + 4 x^3
```

Threads

```mathematica
In[2]:= AlgebraicNumberPolynomial[{2, AlgebraicNumber[Sqrt[2], {1, 2}]}, x]
Out[2]= {2, 1 + 2 x}
```

```mathematica
In[3]:= a = AlgebraicNumber[Sqrt[2 + Sqrt[3]], {1, 2, 3, 4}]; RootReduce[(AlgebraicNumberPolynomial[a, x] /. x -> Sqrt[2 + Sqrt[3]]) == a]
Out[3]= True
```

### Worked examples (2)

```mathematica
In[4]:= AlgebraicNumberPolynomial[2, x]
Out[4]= 2

In[5]:= AlgebraicNumberPolynomial[1/2, x]
Out[5]= 1/2
```

## Algorithm

algebraicnumberpolynomial.c — AlgebraicNumberPolynomial[a, x].

See algebraicnumberpolynomial.h for the contract. The coefficient vector of an AlgebraicNumber object is stored in the object itself, so producing the defining polynomial c0 + c1 x + ... + cn x^n is a purely structural read + build: no field arithmetic and no FLINT. The returned Plus/Times/Power tree is canonicalised by the evaluator on its next fixed-point step (zero terms drop, Times[1, x] folds to x, monomials sort by degree).

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [AlgebraicNumber](../../algebra/AlgebraicNumber/), [RootReduce](../../algebra/RootReduce/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
- Tests: [`tests/test_algebraicnumberpolynomial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumberpolynomial.c)
