# Grad

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Grad[f, {x1, ..., xn}]`**

gives the gradient (D\[f,x1\], ..., D\[f,xn\]) of the scalar f; for an array f a new innermost slot is appended, so a vector field yields its Jacobian.  Equivalent to D\[f, {{x1, ..., xn}}\].

**`Grad[f, {x1, ..., xn}, chart]`**

gives the gradient of a scalar in the orthonormal basis of chart, one of "Cartesian", "Polar", "Cylindrical", "Spherical".

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

```mathematica
In[1]:= Grad[Sin[x^2 + y^2], {x, y}]
Out[1]= {2 x Cos[x^2 + y^2], 2 y Cos[x^2 + y^2]}

In[2]:= Grad[{x y, y z, z x}, {x, y, z}]
Out[2]= {{y, x, 0}, {0, z, y}, {z, 0, x}}

In[3]:= Div[{x^2, y^2, z^2}, {x, y, z}]
Out[3]= 2 x + 2 y + 2 z

In[4]:= Curl[{y, -x}, {x, y}]
Out[4]= -2

In[5]:= Laplacian[x^2 + y^2 + z^2, {x, y, z}]
Out[5]= 6

In[6]:= Div[{r Sin[t], -r Cos[t]}, {r, t}, "Polar"]
Out[6]= 3 Sin[t]

In[7]:= -Grad[k q/r, {r, t, p}, "Spherical"]
Out[7]= {(k q)/r^2, 0, 0}

In[8]:= Laplacian[Sin[r^2], {r, t}, "Polar"] // Simplify
Out[8]= 4 - 4 r^2
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Div](../../calculus/Div/), [Curl](../../calculus/Curl/), [Laplacian](../../calculus/Laplacian/), [D](../../calculus/D/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_vectoranal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_vectoranal.c)
