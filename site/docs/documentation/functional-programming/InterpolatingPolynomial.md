# InterpolatingPolynomial

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`InterpolatingPolynomial[{f1, f2, ...}, x]`**

gives the single polynomial in x reproducing the values fi at x = 1, 2, ..., in nested (Horner) form. With n values the degree is n-1.

**`InterpolatingPolynomial[{{x1, f1}, {x2, f2}, ...}, x]`**

interpolates the values fi at the abscissae xi (arbitrary real, complex, or -- in 1-D -- symbolic).

**`InterpolatingPolynomial[{{{x1, y1, ...}, f1}, ...}, {x, y, ...}]`**

gives the multidimensional interpolating polynomial of lowest total degree.

**`InterpolatingPolynomial[{{xi, fi, dfi, ...}, ...}, x]`**

reproduces derivatives as well as values (the n-th derivative in m-D is a tensor shaped like D\[f, {{x, ...}, n}\]).

<details>
<summary>Notes</summary>

A value or derivative given as Automatic is filled in from the other conditions. The option Modulus -\> n finds the polynomial over the integers modulo n. Exact data give an exact polynomial.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= InterpolatingPolynomial[{1, 4, 9, 16}, x]
Out[1]= 1 + (-1 + x) (1 + x)

In[2]:= InterpolatingPolynomial[{4, 7, 2, {8, 0}, 9}, x]   (* value 8, slope 0 at x=4 *)
Out[2]= 4 + (-1 + x) (3 + (-2 + x) (-4 + (-3 + x) (19/6 + (-4 + x) (-107/36 + 109/72 (-4 + x)))))

In[3]:= Expand[InterpolatingPolynomial[ {{{0, 0}, 1}, {{1, 0}, 7}, {{0, 1}, 10}, {{2, 1}, 40}, {{3, 3}, 151}, {{1, 2}, 47}}, {x, y}]]
Out[3]= 1 + 2 x + 4 x^2 + 3 y + 5 x y + 6 y^2

In[4]:= Expand[InterpolatingPolynomial[ {{-1, Automatic, 0}, {0, 1, 1}, {1, Automatic, 0}}, x]]
Out[4]= 1 + x - 1/3 x^3
```

### Options (1)

```mathematica
In[5]:= InterpolatingPolynomial[{1, 4, 9, 16}, x, Modulus -> 7]
Out[5]= x^2
```

## Algorithm

interp.c

InterpolatingFunction --- piecewise-polynomial interpolation of tabulated data on a regular (tensor-product) grid, plus the Interpolation[] builder. Modelled on Mathematica's InterpolatingFunction object.

```text
  InterpolatingFunction[domain, table]
  InterpolatingFunction[domain, table, ders]
  InterpolatingFunction[domain, table, ders, orders]
  InterpolatingFunction[domain, table, ders, orders, method]

    domain = {{x1min, x1max}, ...}   -- one interval per dimension; the
             number of intervals m is the dimensionality.
    table  = {{coord, val}, ...}                     -- value-only data, or
             {{coord, val, grad, hess, ...}, ...}     -- derivative-supplied.
             coord is a scalar (1-D value-only) or an {x1,...,xm} list.
             grad = D[f,{vars,1}] (length-m vector), hess = D[f,{vars,2}]
             (m x m matrix), etc.
    ders   = {d1, ..., dm}   -- (optional) derivative-of-interpolant orders.
    orders = {o1, ..., om}   -- (optional) interpolation order per dimension.
    method = "Spline" | "Hermite"   -- (optional) interpolation method.
```

Methods (all evaluate the ders-th mixed derivative so D[ifun[..],..] composes):

```text
  default  : sliding-window Newton divided-difference (order min(3,n-1) or the
             requested InterpolationOrder), per dimension, tensor product.
  "Spline" : natural cubic spline (C2; second derivative 0 at the ends),
             tensor product over the full grid.
  "Hermite": tensor-product piecewise cubic Hermite with node slopes estimated
             by 3-point finite differences.
  supplied : derivative-annotated data is interpolated by tensor-product
             Hermite of per-dimension order k = max(K,1) where K is the highest
             supplied derivative order.  Mixed partials that are not supplied
             are filled by central finite differences across the grid.
```

Precision: machine (double) by default; if the data/argument carry MPFR arbitrary precision the MPFR kernels (interp_mpfr.c) are used instead and an EXPR_MPFR is returned.

Builtin ownership: interp_apply / the Interpolation builtin return a fresh Expr* (or NULL to stay unevaluated); inputs are borrowed.

## Implementation notes

**Attributes:** `Protected`.

## See also

[Interpolation](../../functional-programming/Interpolation/), [InterpolatingFunction](../../functional-programming/InterpolatingFunction/), [N](../../arithmetic/N/), [NDArray](../../linear-algebra/NDArray/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_interp_poly.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interp_poly.c)
