### Worked examples

The pseudo-remainder chain starting from the two input polynomials; here it
terminates as soon as a divisor is reached:

```mathematica
In[1]:= SubresultantPolynomialRemainders[x^4 - 1, x^2 - 1, x]
Out[1]= {-1 + x^4, -1 + x^2}
```

For coprime inputs the chain runs all the way down to a nonzero constant:

```mathematica
In[1]:= SubresultantPolynomialRemainders[x^3 - 2 x + 5, x^2 - 3, x]
Out[1]= {5 - 2 x + x^3, -3 + x^2, 5 + x, 22}
```

### Notes

`SubresultantPolynomialRemainders[a, b, x]` returns the remainder chain
`{a, b, R_2, R_3, ...}` obtained by iterating the pseudo-remainder over
`K(coeffs)[x]` until a constant or zero remainder is reached. The final nonzero
entry is the resultant up to content (a constant precisely when `a` and `b` are
coprime). The chain is correct modulo content scaling, which downstream
consumers strip with `primitive[]`; it is the workhorse of the Lazard–Rioboo–
Trager rational-integration pipeline.
