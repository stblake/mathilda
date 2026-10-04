### Worked examples

```mathematica
In[1]:= Exponent[1 + x^2 + x^5, x]  (* the highest power of x *)
```

```mathematica
In[2]:= Exponent[(1 + x)^3, x]  (* expr need not be expanded first *)
```

```mathematica
In[3]:= Exponent[1 + x + x^2, x, List]  (* collect the exponent set with a custom h *)
```

```mathematica
In[4]:= Exponent[a x^2 + b x y^3, {x, y}]  (* a list of forms threads *)
```

```mathematica
In[5]:= Exponent[0, x]  (* the zero polynomial: Max of an empty set *)
```

### Notes

`Exponent[expr, form]` gives the maximum power of `form` in the expanded form of
`expr`, applying `h` (default `Max`) to the set of exponents found;
`Exponent[expr, form, h]` substitutes any other `h`, so `h = List` returns the
whole sorted exponent set. `expr` is expanded first, so it need not be given
expanded, and `form` may be a symbol, a kernel, or a product of terms. The
reading is **purely syntactic** — there is no zero-coefficient recognition — and
the genuine zero polynomial has an empty exponent set, so `Exponent[0, x]` is
`Max[]` = `-Infinity`. Because `Exponent` is `Listable`, a list of forms threads
into a list of exponents.
