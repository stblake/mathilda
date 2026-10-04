### Worked examples

```mathematica
In[1]:= BarnesG[5]  (* the superfactorial 1! 2! 3! = 12 *)
```

```mathematica
In[1]:= BarnesG[Range[6]]  (* threads over the list: G(n) at the first six integers *)
```

```mathematica
In[1]:= BarnesG[9]/BarnesG[8]  (* the functional equation gives G(n+1)/G(n) = Gamma[n] = 7! *)
```

```mathematica
In[1]:= N[BarnesG[7/2], 20]  (* evaluated at a half-integer through the asymptotic continuation *)
```

```mathematica
In[1]:= BarnesG[-2]  (* non-positive integers are (double) zeros *)
```

### Notes

The Barnes G-function satisfies `G(1) = G(2) = 1` and the functional equation
`G(z+1) = Gamma[z] G(z)`. On the positive integers it is the **superfactorial**
`G(n) = prod_{k=1}^{n-2} k!`, computed exactly with GMP; the non-positive
integers are its zeros.

Off the integers there is no elementary closed form, so a non-integer argument
has a value only under `N`, where an asymptotic continuation (twelve fixed
Bernoulli terms plus the Glaisher–Kinkelin constant) carries it to roughly forty
digits. `BarnesG` has no `NDArray` kernel and does not lower under `Compile[]`:
it is exact-integer GMP work or an assembled high-precision continuation, not an
element-wise machine primitive.

`Product` uses `BarnesG` to recognise `prod_{k=1}^{n-1} Gamma[k] = BarnesG[n]`.
