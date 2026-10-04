### Worked examples

```mathematica
In[1]:= BitLength[255]  (* 255 is 11111111 in binary, eight bits *)
```

```mathematica
In[1]:= BitLength[256]  (* one more bit than 255 *)
```

```mathematica
In[1]:= BitLength[2^100]  (* exact for arbitrarily large integers, with no floating point *)
```

```mathematica
In[1]:= BitLength[-256]  (* a negative n uses two's complement: BitLength[BitNot[n]] *)
```

```mathematica
In[1]:= BitLength[{0, 1, 2, 7, 8, 255, 256}]  (* Listable: it threads over a list *)
```

### Notes

`BitLength[n]` is the number of binary bits needed to represent the integer `n`. For `n > 0`
it is an exact `Floor[Log[2, n]] + 1` that never converts through floating point (it uses
GMP's base-2 `mpz_sizeinbase`), so it is exact for arbitrarily large `n`. `BitLength[0]` is
`0`.

For `n < 0` it equals `BitLength[BitNot[n]]`, and in two's complement `BitNot[n] = -n - 1`, so
`BitLength[-1]` is `0`, `BitLength[-2]` is `1`, and `BitLength[-2^k]` is `k`.

`BitLength` is `Protected` and `Listable`, so it threads element-wise over a list. It has a
packed/`NDArray` `int64` fast path and a `Compile[]` lowering, both covering the full `int64`
range including `INT64_MIN` (`BitLength[-2^63]` is `63`). A non-integer argument emits
`BitLength::int` and stays unevaluated; a symbolic argument is left unevaluated silently.
