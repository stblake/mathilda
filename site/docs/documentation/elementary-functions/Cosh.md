# Cosh

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Cosh[z]
    gives the hyperbolic cosine of z, (Exp[z] + Exp[-z]) / 2.
Cosh is Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_cosh` (`src/hyperbolic.c`) mirrors the trig cascade. (1) `strip_inverse_call` folds `Cosh[ArcCosh[x]] -> x`. (2) `try_simp_forward_of_inverse_hyp` handles `Cosh[ArcSinh[x]] -> Sqrt[1+x^2]` and `Cosh[ArcTanh[x]] -> 1/Sqrt[1-x^2]`. (3) `even_fold` for evenness `Cosh[-x] -> Cosh[x]`. (4) `hyp_i_fold` rewrites `Cosh[I y] -> Cos[y]`. (5) `Cosh[0] -> 1`; `Cosh[±Infinity] -> Infinity`. (6) Numeric fallback: MPFR via `mpfr_cosh`/`mpfr_complex_cosh`, else `get_approx` + C99 `ccosh` for inexact real/complex inputs. Otherwise `NULL`. There is no exact rational-multiple-of-Pi table for the hyperbolic heads.

**Data structures.** `Expr*` trees built with the shared `make_*` helpers.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/hyperbolic.c`](https://github.com/stblake/mathilda/blob/main/src/hyperbolic.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Cosh[0]
Out[1]= 1

In[2]:= Cosh[-x]
Out[2]= Cosh[x]
```

```mathematica
In[1]:= Cosh[Pi I]
Out[1]= -1
```

```mathematica
In[1]:= Cosh[ArcSinh[x]]
Out[1]= Sqrt[1 + x^2]
```

```mathematica
In[1]:= TrigExpand[Cosh[a + b]]
Out[1]= Cosh[a] Cosh[b] + Sinh[a] Sinh[b]
```

```mathematica
In[1]:= Series[Cosh[x], {x, 0, 6}]
Out[1]= 1 + 1/2 x^2 + 1/24 x^4 + 1/720 x^6 + O[x]^7

In[2]:= N[Cosh[1], 40]
Out[2]= 1.5430806348152437784779056207570616826015
```

### Notes

`Cosh[z]` is the hyperbolic cosine, `(Exp[z] + Exp[-z])/2`. It is even, so the sign of the argument is dropped. `Cosh` is Listable.
