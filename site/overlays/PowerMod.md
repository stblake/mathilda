---
status: Stable
references:
  - "Knuth, \"The Art of Computer Programming, Vol. 2: Seminumerical Algorithms\", on modular exponentiation and modular inverses."
  - "von zur Gathen & Gerhard, \"Modern Computer Algebra\", on the extended Euclidean algorithm."
---
### Worked examples

```mathematica
In[1]:= PowerMod[3, 1000000, 7]
Out[1]= 4
```

```mathematica
In[1]:= PowerMod[2, 100, 101]
Out[1]= 1
```

```mathematica
In[1]:= PowerMod[7, -1, 11]
Out[1]= 8
```

```mathematica
In[1]:= PowerMod[2, -1, 4]
Out[1]= PowerMod[2, -1, 4]
```

### Notes

`PowerMod[a, b, m]` computes `a^b mod m` using binary modular exponentiation, so a
huge exponent like `PowerMod[3, 1000000, 7] = 4` is evaluated without forming the
full power. `PowerMod[2, 100, 101] = 1` illustrates Fermat's little theorem, since
`101` is prime and `2^100 ≡ 1 (mod 101)`. A negative exponent computes the modular
inverse via the extended Euclidean algorithm: `PowerMod[7, -1, 11] = 8` because
`7*8 = 56 ≡ 1 (mod 11)`. When the inverse does not exist (the base shares a factor
with the modulus, as with `2` and `4`), the call returns unevaluated rather than
producing a spurious result.
