# Divisible

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Divisible[n, m] yields True if n is an integer multiple of m -- the divisibility relation m | n, effectively Mod[n, m] == 0 -- and False otherwise.`**

<details>
<summary>Notes</summary>

Works for machine and BigInt integers, Gaussian integers, rationals, and exact numeric quantities (the quotient n/m must reduce to an integer or Gaussian integer).  Returns False unless n and m are manifestly divisible; symbolic, non-numeric arguments are left unevaluated.  Listable.

</details>

## Examples (18)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

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

### Worked examples (9)

```mathematica
In[7]:= Divisible[10^3000 + 1, 16001]
Out[7]= True

In[8]:= Divisible[0, 0]
Out[8]= True

In[9]:= Divisible[6, 0]
Out[9]= False

In[10]:= Divisible[10, -2]
Out[10]= True

In[11]:= Divisible[3 + I, 1 - I]
Out[11]= True

In[12]:= Divisible[3/2, 1/2]
Out[12]= True

In[13]:= Divisible[2 Pi, Pi/2]
Out[13]= True

In[14]:= Divisible[Sqrt[6], Sqrt[2]]
Out[14]= False

In[15]:= Divisible[{1, 2, 3, 4, 5, 6}, 2]
Out[15]= {False, True, False, True, False, True}
```

### Applications (3)

```mathematica
In[16]:= Divisible[100, 4]
Out[16]= True

In[17]:= Divisible[100, 7]
Out[17]= False

In[18]:= Divisible[10 + 5 I, 1 + 2 I]
Out[18]= True
```

## Implementation notes

**Algorithm.** `builtin_divisible` requires exactly two arguments (otherwise it emits
`Divisible::argm`/`argt` and leaves the call unevaluated). When both `n` and `m` are
integer-like it answers directly with GMP's `mpz_divisible_p` — exact at any precision,
and honouring the convention that divisibility by `0` holds iff `n == 0`. Otherwise it
assembles the quotient `Times[n, Power[m, -1]]`, evaluates it with `eval_and_free`, and
returns `True` iff the result is an integer or a Gaussian integer (so `3 + I` is divisible
by `1 - I`, `3/2` by `1/2`, and `2 Pi` by `Pi/2`). If the quotient is non-integral and
both arguments are concrete numeric quantities it returns `False`; if either argument is
symbolic/non-numeric it returns `NULL`, leaving the call unevaluated so user rules and
pattern matching can apply.

**Data structures.** The integer path uses two `mpz_t`. The general path builds the
quotient as `Power`/`Times` `Expr` trees and runs them back through the evaluator;
`divisible_is_numeric_quantity` walks the result to decide numeric-vs-symbolic, recognising
exact numbers, the named constants (`Pi`, `E`, `EulerGamma`, …), `Complex`/`Rational`, and
any `NumericFunction` applied to numeric quantities. There is no ND/packed/`Compile` path —
`Divisible` returns a Boolean and is `Listable`, threaded by the evaluator before the
builtin runs.

**Complexity / limits.** The integer test is a single quasi-linear `mpz_divisible_p`, so
`Divisible[10^3000 + 1, 16001]` is exact and cheap; the general path costs one full
evaluation of the quotient. Sign is ignored through the multiple test, `Divisible[0, 0]`
is `True`, and a numeric but non-divisible pair is `False`, while symbolic arguments stay
unevaluated.

- Machine integers and GMP bigints: tested directly with `mpz_divisible_p`, so large cases such as `Divisible[10^3000 + 1, 16001]` → `True` are exact. By the GMP convention, divisibility by `0` holds iff `n == 0` (`Divisible[0, 0]` → `True`, `Divisible[6, 0]` → `False`); sign is ignored (`Divisible[10, -2]` → `True`).
- Gaussian integers, rationals, and exact numeric quantities: the quotient `n/m` is formed and evaluated; the result is `True` iff it reduces to an integer or a Gaussian integer. So `Divisible[3 + I, 1 - I]` → `True`, `Divisible[3/2, 1/2]` → `True`, `Divisible[2 Pi, Pi/2]` → `True`, while `Divisible[Sqrt[6], Sqrt[2]]` → `False`.
- `Listable`: threads element-wise over lists, e.g. `Divisible[{1, 2, 3, 4, 5, 6}, 2]` → `{False, True, False, True, False, True}`.
- Symbolic, non-numeric arguments leave the call unevaluated (e.g. `Divisible[x, 2]`).
- Diagnostics: too few arguments emit `Divisible::argm`, too many emit `Divisible::argt`; both leave the call unevaluated.

**Attributes:** `Listable`, `Protected`.

## References

- G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — divisibility (Chapter I).
- Source: [`src/numbertheory/divisible.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory/divisible.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_divisible.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisible.c)
- Tests: [`tests/test_divisors.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisors.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### The divisibility relation

`Divisible[n, m]` tests the relation `m ∣ n` — whether `n` is an integer multiple of `m`,
equivalently whether `Mod[n, m] == 0`. It extends beyond the ordinary integers: over the
Gaussian integers `Z[i]`, `m ∣ n` when the quotient `n/m` is itself a Gaussian integer.
