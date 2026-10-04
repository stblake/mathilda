# Log10

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Log10[z]`**

gives the base-10 logarithm of z, Log\[10, z\] = Log\[z\] / Log\[10\].

<details>
<summary>Notes</summary>

Exact powers of 10 give exact results (Log10\[1000\] = 3); a symbolic z gives Log\[z\]/Log\[10\]. Listable; maps packed arrays and compiles.

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

Powers of ten come back exact, via Log[10, z]'s power detection

```mathematica
In[7]:= Log10[1000]
Out[7]= 3
```

Listable, so a whole list maps at once

```mathematica
In[8]:= Log10[{1, 10, 100, 1000}]
Out[8]= {0, 1, 2, 3}
```

An arbitrary-precision value routes through the MPFR logarithm

```mathematica
In[9]:= N[Log10[2], 30]
Out[9]= 0.3010299956639811952137388947246
```

A negative real lands off the real axis, exactly as Log does

```mathematica
In[10]:= Log10[-1.]
Out[10]= 0.0 + 1.36438*I
```

Differentiated through the Log[10, z] = Log[z]/Log[10] definition

```mathematica
In[11]:= D[Log10[x], x]
Out[11]= 1/(Log[10] x)
```

## Implementation notes

**Algorithm.** `builtin_log10` is `fixed_base_log(res, "Log10", 10, log10)`.
Log10[z] is defined, as in Mathematica, as Log[10, z] = Log[z]/Log[10]. One fast
arm handles a positive finite machine real directly with libm's `log10`, because
the quotient of two separately rounded logarithms is wrong by an ulp at exact
powers — `log(1000.)/log(10.)` is 2.9999999999999996 where `log10(1000.)` is
exactly 3. Every other argument — an exact power of ten (detected exactly), a
symbolic z, a negative or complex argument, an arbitrary-precision number, or
zero — is forwarded to `Log[10, z]`, so the two spellings can never disagree.

**Data structures.** A single `double` on the fast arm; otherwise a two-argument
`Log[10, z]` expression handed back to the evaluator. The ND kernel
(`UK_FIXED_LOG(Log10, log10, ND_LN10)`, registered with `REG_U`) maps a packed or
visible real `NDArray` element-wise through libm `log10`; a negative or complex
element escapes to `clog(z)/ln(10)`, so exact powers of ten stay exact on the
buffer as well.

**Complexity / limits.** `O(1)` per element; `Compile[]` lowers it at scalar and
rank-1 shapes. `Log10[0]` is `-Infinity` and the negative real axis returns the
principal complex value, both inherited from `Log`.

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

**See also:** [Log2](../../elementary-functions/Log2/)

- Source: [`src/logexp.c`](https://github.com/stblake/mathilda/blob/main/src/logexp.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_logexp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_logexp.c)
- Tests: [`tests/test_numeric.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric.c)

## Notes & additional examples

### Notes

`Log10[z]` is `Log[10, z] = Log[z]/Log[10]`, so every rule `Log` knows — exact
powers, the branch cut on the negative real axis, arbitrary precision — is
inherited. The only thing `Log10` adds is the dedicated libm `log10` on a
positive machine real, which keeps `Log10[1000.]` exactly `3.` where the quotient
of two rounded logarithms would land an ulp away.

The real-valued `NDArray` kernel maps a packed buffer element-wise and compiles
at scalar and rank-1 shapes, so `Log10` in a `Compile[]`d body or over a packed
column runs on the buffer rather than one boxed number at a time.
