### Worked examples

A trivial value — Euler's solution of the Basel problem, ζ(2) = π²/6:

```mathematica
In[1]:= Zeta[2]
Out[1]= 1/6 Pi^2
```

Even positive integers are rational multiples of powers of π (via Bernoulli numbers):

```mathematica
In[1]:= Zeta[6]
Out[1]= 1/945 Pi^6
```

Negative-integer values are rational, exposing the Bernoulli-number connection ζ(−n) = −Bₙ₊₁/(n+1); the negative even integers are the trivial zeros:

```mathematica
In[1]:= Zeta[-1]
Out[1]= -1/12

In[1]:= Zeta[-3]
Out[1]= 1/120

In[1]:= Table[Zeta[-2 n], {n, 1, 4}]
Out[1]= {0, 0, 0, 0}
```

Apéry's constant ζ(3) to 40 digits, evaluated through the MPFR kernel:

```mathematica
In[1]:= N[Zeta[3], 40]
Out[1]= 1.2020569031595942853997381615114499907651
```

A value on the critical line at the first nontrivial zero ρ ≈ 1/2 + 14.134725 i — the result is numerically zero to within the input precision, confirming the zero of the Riemann hypothesis:

```mathematica
In[1]:= N[Zeta[1/2 + 14.134725 I], 10]
Out[1]= 1.7674298414e-08 - 1.1102028931e-07*I
```

The Laurent expansion about the pole s = 1 defines the Stieltjes constants, with residue 1 and constant term the Euler–Mascheroni constant γ:

```mathematica
In[1]:= Series[Zeta[s], {s, 1, 2}]
Out[1]= 1/(s - 1) + EulerGamma + -StieltjesGamma[1] (s - 1) + 1/2 StieltjesGamma[2] (s - 1)^2 + O[s - 1]^3
```

The Hurwitz zeta at an integer second argument reduces to ζ(s) minus a finite power sum, mixing exact π-power and rational parts:

```mathematica
In[1]:= Zeta[4, 5]
Out[1]= -22369/20736 + 1/90 Pi^4
```
