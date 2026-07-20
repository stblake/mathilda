# Fourier

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Fourier[list]
    gives the discrete Fourier transform of a list of complex numbers.
Fourier[list, {p1, p2, ...}]
    returns the specified positions of the discrete Fourier transform.

The transform of a length-n list u is v[s] = 1/n^((1-a)/2) Sum_r u[r]
Exp[2 Pi I b (r-1)(s-1)/n], with {a,b} set by the FourierParameters
option (default {0,1}; {-1,1} data analysis, {1,-1} signal processing).
Exact input is first numericalised with N; the list may be a nested
rectangular array for a multidimensional transform. Symbolic input
yields the exact transform in terms of roots of unity.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Fourier[{1, 1, 2, 2, 1, 1, 0, 0}]
Out[1]= {2.82843, -0.5 + 1.20711*I, 0.0, 0.5 - 0.207107*I, 0.0, 0.5 + 0.207107*I, 0.0, -0.5 - 1.20711*I}

In[2]:= Abs[Fourier[{1, 2, 3, 4, 5, 6}]]^2
Out[2]= {73.5, 6.0, 2.0, 1.5, 2.0, 6.0}

In[3]:= Fourier[{a, b, c, d}]
Out[3]= {1/2 (a + b + c + d), 1/2 (a + I b - c - I d), 1/2 (a - b + c - d), 1/2 (a - I b - c + I d)}

In[4]:= Fourier[{1, 0, 1, 0, 0, 1, 0, 0, 0, 1}, FourierParameters -> {-1, 1}][[1]]
Out[4]= 0.4
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/fourier-transforms.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/fourier-transforms.md)
