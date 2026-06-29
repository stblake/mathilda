# Divisible

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Divisible[n, m]
    yields True if n is divisible by m, and False otherwise.
n is divisible by m when n is an integer multiple of m; this is
effectively Mod[n, m] == 0.  Works for machine and BigInt integers,
Gaussian integers, rationals, and exact numeric quantities (the
quotient n/m must reduce to an integer or Gaussian integer).  Returns
False unless n and m are manifestly divisible; symbolic, non-numeric
arguments are left unevaluated.  Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Divisible[10, 2]
Out[1]= True

In[2]:= Divisible[5, 2]
Out[2]= False

In[3]:= Divisible[3 + I, 1 - I]
Out[3]= True

In[4]:= Divisible[2 Pi, Pi/2]
Out[4]= True

In[5]:= Divisible[Sqrt[6], Sqrt[2]]
Out[5]= False

In[6]:= Divisible[{1, 2, 3, 4, 5, 6}, 2]
Out[6]= {False, True, False, True, False, True}
```

## Implementation notes

- Machine integers and GMP bigints: tested directly with `mpz_divisible_p`, so large cases such as `Divisible[10^3000 + 1, 16001]` → `True` are exact. By the GMP convention, divisibility by `0` holds iff `n == 0` (`Divisible[0, 0]` → `True`, `Divisible[6, 0]` → `False`); sign is ignored (`Divisible[10, -2]` → `True`).
- Gaussian integers, rationals, and exact numeric quantities: the quotient `n/m` is formed and evaluated; the result is `True` iff it reduces to an integer or a Gaussian integer. So `Divisible[3 + I, 1 - I]` → `True`, `Divisible[3/2, 1/2]` → `True`, `Divisible[2 Pi, Pi/2]` → `True`, while `Divisible[Sqrt[6], Sqrt[2]]` → `False`.
- `Listable`: threads element-wise over lists, e.g. `Divisible[{1, 2, 3, 4, 5, 6}, 2]` → `{False, True, False, True, False, True}`.
- Symbolic, non-numeric arguments leave the call unevaluated (e.g. `Divisible[x, 2]`).
- Diagnostics: too few arguments emit `Divisible::argm`, too many emit `Divisible::argt`; both leave the call unevaluated.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
