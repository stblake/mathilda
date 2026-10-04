### Worked examples

```mathematica
In[1]:= DigitSum[12345]  (* sum of the base-10 digits *)
```

```mathematica
In[1]:= DigitSum[255, 2]  (* base 2: 255 is 11111111 *)
```

```mathematica
In[1]:= DigitSum[255, 16]  (* base 16: FF is 15 + 15 *)
```

```mathematica
In[1]:= DigitSum[-987]  (* the sign is discarded *)
```

```mathematica
In[1]:= DigitSum[2^100]  (* arbitrary-precision bignum, summed in GMP *)
```

```mathematica
In[1]:= DigitSum[{123, 4567, 89}]  (* Listable: threads over the list *)
```

```mathematica
In[1]:= Mod[n - DigitSum[n], 9] /. n -> 123456  (* casting out nines: n is DigitSum[n] mod 9 *)
```

### Notes

`DigitSum[n]` adds up the base-10 digits of `n`; `DigitSum[n, b]` uses base
`b >= 2`. The sign of `n` is always discarded (`DigitSum[-987]` equals
`DigitSum[987]`), and `DigitSum[0]` is `0`. It is exactly
`Total[IntegerDigits[n, b]]`, computed without building the digit list.

Both machine integers and arbitrary-precision bignums work uniformly, since the
digit extraction runs in GMP. The casting-out-nines identity `n ≡ DigitSum[n]
(mod b-1)` (here `mod 9` for base 10) follows directly from summing digits.

`DigitSum` is `Listable`, so it maps over the elements of a list; it is **not**
accelerated on a packed buffer, and it does not lower in `Compile[]`. Its close
relatives are `IntegerDigits` (the digit list), `DigitCount` (per-value digit
tallies), and `IntegerLength` (how many digits there are).
