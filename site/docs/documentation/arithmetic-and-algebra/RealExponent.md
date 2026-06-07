# RealExponent

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
RealExponent[x] gives Log[10, |x|] -- the base-10 real exponent of x.
RealExponent[x, b] gives Log[b, |x|] in the specified base b.
Accepts Integer, BigInt, Rational, Real, and (with USE_MPFR) MPFR inputs, plus symbolic numeric values such as Pi, E, or Pi^Pi.  Result is a machine Real unless an MPFR input lifts it to MPFR at that precision.  Exact zero gives -Infinity; machine 0. gives Log[b, $MinMachineNumber] (~ -307.65 in base 10); MPFR 0 with precision p digits gives -p / Log10[b].  Threads over lists.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= RealExponent[123.456]
Out[1]= 2.09151

In[2]:= RealExponent[123.456, 2]
Out[2]= 6.94785

In[3]:= RealExponent[N[Pi, 32]]
Out[3]= 0.497149872694133854351268288290899

In[4]:= RealExponent[Pi, E]
Out[4]= 1.14473

In[5]:= RealExponent[987654321/123456789]
Out[5]= 0.90309

In[6]:= RealExponent[{1, 2, 3, 4, 5}]
Out[6]= {0.0, 0.30103, 0.477121, 0.60206, 0.69897}

In[7]:= Table[RealExponent[Pi, b], {b, {2, 3, 5, 7, 10}}]
Out[7]= {1.6515, 1.04198, 0.711261, 0.588275, 0.49715}

In[8]:= RealExponent[0]
Out[8]= -Infinity
```

## Implementation notes

`builtin_real_exponent` returns `RealExponent[x]` / `RealExponent[x, b]` — essentially `⌊Log_b|x|⌋`, the exponent of the leading digit. It rejects true (non-zero-imaginary) `Complex` inputs (`RealExponent::realx`/`::ibase`) and bad arg counts (`RealExponent::argt`). Symbolic constants (Pi, E, …) and either argument are numericalised to a recognised numeric kind at a working precision lifted to cover any MPFR input (`+32` guard bits, so the downstream `Log` keeps precision), then the floor of the base-b logarithm of `|x|` is taken.

- `Protected`, `Listable`. Threads over lists in any argument position.
- Accepts `Integer`, `BigInt`, `Rational`, machine `Real`, and (under

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/real.c`](https://github.com/stblake/mathilda/blob/main/src/real.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
