# Cosh

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Cosh[z]`**

gives the hyperbolic cosine of z, (Exp\[z\] + Exp\[-z\]) / 2.

<details>
<summary>Notes</summary>

Cosh is Listable.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (7)

```mathematica
In[1]:= Cosh[0]
Out[1]= 1

In[2]:= Cosh[-x]
Out[2]= Cosh[x]

In[3]:= Cosh[Pi I]
Out[3]= -1

In[4]:= Cosh[ArcSinh[x]]
Out[4]= Sqrt[1 + x^2]

In[5]:= TrigExpand[Cosh[a + b]]
Out[5]= Cosh[a] Cosh[b] + Sinh[a] Sinh[b]

In[6]:= Series[Cosh[x], {x, 0, 6}]
Out[6]= 1 + 1/2 x^2 + 1/24 x^4 + 1/720 x^6 + O[x]^7

In[7]:= N[Cosh[1], 40]
Out[7]= 1.5430806348152437784779056207570616826015
```

## Implementation notes

**Algorithm.** `builtin_cosh` (`src/hyperbolic.c`) mirrors the trig cascade. (1) `strip_inverse_call` folds `Cosh[ArcCosh[x]] -> x`. (2) `try_simp_forward_of_inverse_hyp` handles `Cosh[ArcSinh[x]] -> Sqrt[1+x^2]` and `Cosh[ArcTanh[x]] -> 1/Sqrt[1-x^2]`. (3) `even_fold` for evenness `Cosh[-x] -> Cosh[x]`. (4) `hyp_i_fold` rewrites `Cosh[I y] -> Cos[y]`. (5) `Cosh[0] -> 1`; `Cosh[±Infinity] -> Infinity`. (6) Numeric fallback: MPFR via `mpfr_cosh`/`mpfr_complex_cosh`, else `get_approx` + C99 `ccosh` for inexact real/complex inputs. Otherwise `NULL`. There is no exact rational-multiple-of-Pi table for the hyperbolic heads.

**Data structures.** `Expr*` trees built with the shared `make_*` helpers.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/hyperbolic.c`](https://github.com/stblake/mathilda/blob/main/src/hyperbolic.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_besseli.c`](https://github.com/stblake/mathilda/blob/main/tests/test_besseli.c)
- Tests: [`tests/test_complexexpand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_complexexpand.c)
- Tests: [`tests/test_coshintegral.c`](https://github.com/stblake/mathilda/blob/main/tests/test_coshintegral.c)
- Tests: [`tests/test_deriv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_deriv.c)

## Notes & additional examples

### Notes

`Cosh[z]` is the hyperbolic cosine, `(Exp[z] + Exp[-z])/2`. It is even, so the sign of the argument is dropped. `Cosh` is Listable.
