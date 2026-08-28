---
references:
  - "G. E. Andrews, *The Theory of Partitions*, Cambridge University Press, 1998 — distinct-parts partitions and Euler's identity."
---

### Distinct parts equal odd parts

`PartitionsQ[n]` counts the partitions of `n` into **distinct** parts. Euler's celebrated
identity says this equals the number of partitions of `n` into **odd** parts — the
generating-function identity `∏ (1 + x^k) = ∏ 1/(1 - x^(2k-1))`. Contrast
[`PartitionsP`](PartitionsP.md), which counts *all* partitions.

### Worked examples

```mathematica
In[1]:= PartitionsQ[10]
Out[1]= 10
```

```mathematica
In[1]:= PartitionsQ[6]
Out[1]= 4
```

The four partitions of `6` into distinct parts — `6`, `5+1`, `4+2`, `3+2+1` — match the four
into odd parts — `5+1`, `3+3`, `3+1+1+1`, `1+1+1+1+1+1` — as Euler's theorem promises.
