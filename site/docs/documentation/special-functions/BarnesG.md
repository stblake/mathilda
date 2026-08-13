# BarnesG

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BarnesG[z]`**

gives the Barnes G-function.

<details>
<summary>Notes</summary>

G(z+1) = Gamma\[z\] G(z) with G(1)=G(2)=1; for a positive integer n, G(n+1) = prod\_{k=1}^{n-1} k! (exact via GMP), and G(m)=0 for non-positive integer m. A non-integer numeric order (under N) evaluates from the Barnes asymptotic expansion plus the Gamma recurrence (real, complex, arbitrary precision); symbolic orders stay unevaluated. Listable, NumericFunction.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= BarnesG[5]
Out[1]= 12

In[2]:= N[BarnesG[6.0]]
Out[2]= 288.0

In[3]:= N[BarnesG[13/2], 30]
Out[3]= 2548.745769568498989735906104648

In[4]:= N[BarnesG[2.5 + 1.0 I]]
Out[4]= 0.743798 - 0.0953168*I
```

### Worked examples (1)

```mathematica
In[5]:= Product[Gamma[i], {i, 1, n-1}]
Out[5]= BarnesG[n]
```

## Algorithm

Mathilda -- BarnesG[z], the Barnes G-function.

```text
  G(1) = G(2) = 1,   G(z+1) = Gamma[z] G(z),
  integer:  G(n+1) = prod_{k=1}^{n-1} k!   (the superfactorial),
            G(m) = 0 for non-positive integer m (double zeros).
```

Exact for integer orders (GMP); non-integer orders are left unevaluated (the

```text
LogGamma/zeta'(-1) asymptotic continuation is not implemented).  N at an
integer order routes through the exact value and numericalize.  Used by
```

Product to recognise prod_{k=1}^{n-1} Gamma[k] = BarnesG[n].

Memory: honours the builtin ownership contract.

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Product](../../calculus/Product/), [N](../../arithmetic/N/), [Gamma](../../special-functions/Gamma/), [Log](../../elementary-functions/Log/), [Exp](../../elementary-functions/Exp/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_ndsolve_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndsolve_compile.c)
