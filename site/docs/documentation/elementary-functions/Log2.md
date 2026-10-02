# Log2

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Log2[z]`**

gives the base-2 logarithm of z, Log\[2, z\] = Log\[z\] / Log\[2\].

<details>
<summary>Notes</summary>

Exact powers of 2 give exact results (Log2\[1024\] = 10); a symbolic z gives Log\[z\]/Log\[2\]. Listable; maps packed arrays and compiles.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= Log10[100]
Out[1]= 2

In[2]:= Log10[2.]
Out[2]= 0.30103

In[3]:= Log10[x]
Out[3]= Log[x]/Log[10]

In[4]:= Log10[1/1000]
Out[4]= -3

In[5]:= Log2[1024]
Out[5]= 10

In[6]:= Log2[{1., 2., 8.}]
Out[6]= {0.0, 1.0, 3.0}
```

## Implementation notes

- `Listable`, `NumericFunction`, `Protected`, matching Mathematica.
- Defined as `Log[10, z]` / `Log[2, z]`, so exact powers of the base are
  exact (`Log10[1000] = 3`, `Log10[1/100] = -2`, `Log2[1/8] = -3`) and a
  symbolic or non-power argument gives Mathematica's `Log[z]/Log[10]` form.
  Zero, the negative axis, complex and arbitrary-precision arguments follow
  `Log[b, z]` (`Log10[0] = -Infinity`, `Log10[0.] = Indeterminate`,
  `Log10[-1.] = 0. + 1.36438 I`).
- A positive machine real is evaluated with libm `log10` / `log2` rather than
  `log(z)/log(b)`, so an exact power stays exact: `Log10[1000.]` is `3.`, not
  `2.9999999999999996`.
- **Fast paths.** Each has an escaping NDArray kernel (libm `log10`/`log2` on
  the positive axis, `clog(z)/Log[b]` elsewhere, promoting to complex only
  when an element leaves the real axis), so `Log10` threads over a visible
  `NDArray[...]` and a packed real array stays packed. Both lower in
  `Compile[]` (hence auto-compile) at scalar and rank-1 array shape; the
  scalar lowering bails to the interpreter off the positive axis.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Log10](../../elementary-functions/Log10/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_logexp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_logexp.c)
