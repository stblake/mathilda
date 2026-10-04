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

## Examples (11)

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

### Applications (5)

Exact powers of two stay exact

```mathematica
In[7]:= Log2[1024]
Out[7]= 10
```

Listable over a vector of powers

```mathematica
In[8]:= Log2[{1, 2, 4, 8, 16}]
Out[8]= {0, 1, 2, 3, 4}
```

The base-2 logarithm of ten to thirty digits

```mathematica
In[9]:= N[Log2[10], 30]
Out[9]= 3.32192809488736234787031942949
```

A positive machine real goes straight to libm log2

```mathematica
In[10]:= Log2[0.1]
Out[10]= -3.32193
```

Derivative through the Log[2, z] definition

```mathematica
In[11]:= D[Log2[x], x]
Out[11]= 1/(Log[2] x)
```

## Implementation notes

**Algorithm.** `builtin_log2` is `fixed_base_log(res, "Log2", 2, log2)`, the exact
base-2 twin of `Log10`. Log2[z] = Log[2, z] = Log[z]/Log[2]. A positive finite
machine real goes straight to libm's `log2` — the one-ulp error of
`log(z)/log(2)` at exact powers of two (so `Log2[1024.]` would miss `10.`) is why
the dedicated libm call is used rather than a ratio of logarithms. Exact powers
of two, symbolic, negative, complex, arbitrary-precision, and zero arguments all
fall through to `Log[2, z]`, keeping the two spellings consistent.

**Data structures.** A `double` on the fast arm, else a two-argument `Log[2, z]`
tree. The ND kernel (`UK_FIXED_LOG(Log2, log2, ND_LN2)`, `REG_U`) runs a packed or
visible real `NDArray` element-wise through libm `log2`, escaping a negative or
complex element to `clog(z)/ln(2)`.

**Complexity / limits.** `O(1)` per element; `Compile[]` lowers it at scalar and
rank-1 shapes. `Log2[0]` is `-Infinity`; the negative axis returns the principal
complex value, both from `Log`.

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

- Source: [`src/logexp.c`](https://github.com/stblake/mathilda/blob/main/src/logexp.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_logexp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_logexp.c)

## Notes & additional examples

### Notes

`Log2[z]` is `Log[2, z] = Log[z]/Log[2]`, the base-2 counterpart of `Log10`.
Information-theoretic and computer-science uses — bits of entropy, tree depth,
complexity exponents — want a base-2 logarithm, and spelling it `Log2` both reads
clearly and avoids the one-ulp drift of `Log[z]/Log[2]` at exact powers of two.

Like `Log10`, it carries an `NDArray` kernel and a `Compile[]` lowering, so a
packed buffer or a compiled body evaluates `Log2` on the buffer directly.
