# BitLength

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BitLength[n] gives the number of binary bits needed to represent the integer n.`**

<details>
<summary>Notes</summary>

For n \> 0, BitLength\[n\] is Floor\[Log\[2, n\]\] + 1; BitLength\[0\] is 0. For n \< 0, BitLength\[n\] is equivalent to BitLength\[BitNot\[n\]\].

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= BitLength[255]
Out[1]= 8

In[2]:= BitLength[256]
Out[2]= 9

In[3]:= BitLength[2^100]
Out[3]= 101

In[4]:= BitLength[100!]
Out[4]= 525

In[5]:= BitLength[{0, 1, 2, 7, 8, 255, 256}]
Out[5]= {0, 1, 2, 3, 4, 8, 9}

In[6]:= BitLength[{-1, -2, -8, -256}]
Out[6]= {0, 1, 3, 8}
```

### Scope (1)

```mathematica
In[7]:= BitLength[1, 2] BitLength::argx: BitLength called with 2 arguments; 1 argument is expected.
```

## Options & behaviour

**Diagnostics**. A non-integer argument emits `BitLength::int` and echoes the
call back unevaluated; a wrong-arity call emits `BitLength::argx`. A symbolic
argument is left unevaluated silently.

## Algorithm

bitlength.c

BitLength[n] -- the number of binary bits needed to represent the integer n.

```text
  n > 0 : Floor[Log[2, n]] + 1   (= mpz_sizeinbase(n, 2)).
  n = 0 : 0.
  n < 0 : BitLength[BitNot[n]], and BitNot[n] = -n - 1 in two's complement,
          so BitLength[-1] = 0, BitLength[-2] = 1, BitLength[-2^k] = k.
```

BitLength is Listable (threads over lists automatically). The packed / NDArray fast path is the int64 kernel ndk_BitLength_ii in src/ndinteger.c and the Compile[] lowering is OP_BLEN_I in src/compile/. All arithmetic here is done in GMP, so machine integers and arbitrary-precision bignums are handled uniformly.

## Implementation notes

- `Protected`, `Listable`. Threads element-wise over a list of integers, e.g.
  `BitLength[{0, 1, 2, 7, 8, 255, 256}]`.
- For `n > 0`, `BitLength[n]` is an efficient exact version of
  `Floor[Log[2, n]] + 1` -- it never converts through floating-point and is
  exact for arbitrarily large `n` (via GMP's `mpz_sizeinbase` in base 2).
- `BitLength[0]` is `0`.
- For `n < 0`, `BitLength[n]` is equivalent to `BitLength[BitNot[n]]`. In two's
  complement `BitNot[n] = -n - 1`, so `BitLength[-1]` is `0`, `BitLength[-2]` is
  `1`, and `BitLength[-2^k]` is `k`.
- Works for both machine integers and bignums. The packed / `NDArray` fast path
  (`src/ndinteger.c`) and the `Compile[]` lowering (`OP_BLEN_I`) both handle the
  full `int64` range, including `INT64_MIN` (`BitLength[-2^63]` is `63`), which
  the complement-based kernel computes without overflow.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [NDArray](../../linear-algebra/NDArray/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/bitwise.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/bitwise.md)
- Tests: [`tests/test_bitwise.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bitwise.c)
