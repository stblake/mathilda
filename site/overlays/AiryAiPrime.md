### Worked examples

```mathematica
In[1]:= AiryAiPrime[0]
Out[1]= -1/(3^(1/3) Gamma[1/3])
```

```mathematica
In[1]:= N[AiryAiPrime[0], 40]
Out[1]= -0.25881940379280679840518356018920396347907
```

```mathematica
In[1]:= D[AiryAiPrime[z], z]
Out[1]= z AiryAi[z]
```

```mathematica
In[1]:= AiryAiPrime[1.0 + 1.0 I]
Out[1]= -0.130628 + 0.163068*I
```

### Notes

`AiryAiPrime[z]` is the derivative `Ai'(z)`. Its exact origin value is
`-1/(3^(1/3) Gamma[1/3])`, and differentiating once more recovers the Airy
equation in the form `D[AiryAiPrime[z], z] == z AiryAi[z]`. Complex arguments
evaluate to machine precision (and to arbitrary precision under `N[..., n]`);
`AiryAiPrime` is Listable.
