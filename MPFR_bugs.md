# MPFR Bugs

Known defects and divergences from Mathematica in the arbitrary-precision
(MPFR) subsystem. Each entry records a minimal reproducer, the expected
Mathematica behaviour, the actual Mathilda behaviour, and what is known about
the root cause. Surfaced by the `*_tests` suites and ad-hoc REPL probes.

---

## 1. Printer suppresses trailing zeros for arbitrary-precision reals

**Status:** open · display-only (no precision loss) · found 2026-06-07

**Reproducers**

```
In[1]:= N[0.5, 50]
Out[1]= 0.5                         (* expected: 0.50000000000000000000000000000000000000000000000000 *)

In[2]:= N[2, 50]
Out[2]= 2.0                         (* expected: 2.0000000000000000000000000000000000000000000000000 *)

In[3]:= Re[Log[Exp[Complex[N[0.5, 50], N[0.25, 50]]]]]
Out[3]= 0.5                         (* expected: 0.5000...  (50 significant digits) *)
```

Surfaced as a hard failure in `numeric_tests`:

```
FAIL: Re[Log[Exp[Complex[N[0.5, 50], N[0.25, 50]]]]]
  Expected prefix: 0.5000000000000000000000000000000000000000000000000
  Actual:          0.5
```

**Diagnosis — this is a printing bug, not precision loss.** The value carries
its full declared precision; only the rendered form is wrong:

```
In[4]:= Precision[N[0.5, 50]]
Out[4]= 50.272                      (* full 50-digit precision retained *)

In[5]:= Precision[0.5]
Out[5]= MachinePrecision            (* the bare literal is machine precision, prints 0.5 correctly *)

In[6]:= Log[Exp[N[0.5, 50]]] - N[0.5, 50]
Out[6]= 0.0                         (* round-trips exactly: no value error *)
```

Non-round mantissas already print every significant digit, confirming the
defect is specific to values whose decimal expansion terminates early:

```
In[7]:= N[Pi, 50]
Out[7]= 3.1415926535897932384626433832795028841971693993751   (* correct *)

In[8]:= N[1/3, 50]
Out[8]= 0.333333333333333333333333333333333333333333333333332  (* correct *)
```

**Root cause (suspected):** the `EXPR_MPFR` real formatting path in
`src/print.c` trims trailing zeros from the mantissa. For machine reals that
is correct Mathematica behaviour, but for an arbitrary-precision number the
trailing zeros are *significant* — they communicate the precision — and must
be retained (padded to the significant-digit count derived from `Precision`).
The fix belongs in the printer alone; the numeric kernel is sound. Because
`Re`/`Im` of a `Complex[MPFR, MPFR]` route component output through the same
real formatter, the complex case (`Out[3]`) is the same bug.
