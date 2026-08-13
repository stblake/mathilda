# MoebiusMu

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MoebiusMu[n] gives the Moebius function mu(n): 0 if n has a squared prime factor, otherwise (-1)^k where k is the number of distinct primes. A non-real Gaussian-integer argument is handled over Z[i].`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= MoebiusMu[11]
Out[1]= -1

In[2]:= MoebiusMu[10]
Out[2]= 1

In[3]:= MoebiusMu[1440]
Out[3]= 0

In[4]:= MoebiusMu[{4, 10, 17, 20}]
Out[4]= {0, 1, -1, 0}

In[5]:= MoebiusMu[10^50 + 1]
Out[5]= -1

In[6]:= MoebiusMu[5 + 6 I]
Out[6]= -1
```

## Options & behaviour

> **Packed arrays.** Runs on an `int64` buffer. `MoebiusMu[0]` is undefined,
> so an array containing `0` takes the ordinary path and leaves that element
> unevaluated exactly as the scalar does.

## Implementation notes

- `Listable`, `Protected`.
- Computed directly from the prime factorisation (machine integers and GMP
  bigints handled uniformly); the result is always `0`, `1`, or `-1`.
- The sign of `n` is ignored (`mu(-n) = mu(n)`).
- A non-real Gaussian-integer argument `Complex[a, b]` is handled over `Z[i]`:
  the input is factored into Gaussian primes (the unit factor does not count),
  giving `0` for a repeated Gaussian prime factor and `(-1)^m` otherwise.
- Non-integer or zero `n` is left unevaluated; a wrong argument count issues a
  `MoebiusMu::argx` message.

**Attributes:** `Listable`, `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_moebiusmu.c`](https://github.com/stblake/mathilda/blob/main/tests/test_moebiusmu.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_primenu.c`](https://github.com/stblake/mathilda/blob/main/tests/test_primenu.c)
