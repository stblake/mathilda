# NumberForm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NumberForm[expr, n] prints approximate real numbers in expr to n-digit`**

precision; NumberForm\[expr, {n, f}\] uses n digits with f to the right of the decimal point; NumberForm\[expr\] uses the default options. Works on integers as well.  It is an inert print wrapper: the head remains in the expression and only changes how it is displayed. Options: DefaultPrintPrecision, DigitBlock, ExponentFunction, ExponentStep, NumberFormat, NumberMultiplier, NumberPadding, NumberPoint, NumberSeparator, NumberSigns, ScientificNotationThreshold, SignPadding.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= NumberForm[N[Pi], 10]
Out[1]= 3.141592654
```

### Options (2)

```mathematica
In[2]:= NumberForm[10^9, DigitBlock -> 3]
Out[2]= 1,000,000,000

In[3]:= NumberForm[{8.^5, 11.^7, 13.^9}, NumberFormat -> (Row[{#1, "e", #3}] &)]
Out[3]= {32768.e, 1.94872e7, 1.06045e10}
```

### Applications (5)

Ten significant digits

```mathematica
In[4]:= NumberForm[N[Pi], 10]
Out[4]= 3.141592654
```

Thousands separators

```mathematica
In[5]:= NumberForm[10^9, DigitBlock -> 3]
Out[5]= 1,000,000,000
```

Five significant figures, two decimals

```mathematica
In[6]:= NumberForm[1234.567, {5, 2}]
Out[6]= 1234.60
```

Custom layout

```mathematica
In[7]:= NumberForm[{8.^5, 11.^7, 13.^9}, NumberFormat -> (Row[{#1, "e", #3}] &)]
Out[7]= {32768.e, 1.94872e7, 1.06045e10}
```

The wrapper head survives in the tree

```mathematica
In[8]:= FullForm[NumberForm[1.23, 2]]
Out[8]= NumberForm[1.23, 2]
```

## Algorithm

```text
 Mathilda — NumberForm and a minimal Row.  See numberform.h for the design.
```

NumberForm is a PRINT WRAPPER: builtin_numberform is inert (returns NULL) so the head survives in the tree, and all display work happens here, driven by print.c which installs an active NumberFormCtx and routes every numeric leaf through numberform_render_number.

The per-number pipeline (nf_format_parts):

```text
  1. reject non-finite; special-case zero and exact integers.
  2. extract the value's `count` significant base-10 digits + decimal
     exponent (mpfr_get_str for MPFR, "%.*e" for machine reals).
  3. decide scientific vs decimal (ScientificNotationThreshold, or the
     caller's ExponentFunction), and the displayed exponent (ExponentStep).
  4. place the digits into integer/fractional strings; {n,f} re-rounds to f
     fractional digits, plain-n drops trailing zeros.
  5. apply DigitBlock grouping, then assemble sign / point / multiplier /
     NumberFormat.
```

The same nf_format_parts drives the measure pass, so alignment widths and the printed output can never disagree.

## Implementation notes

**Algorithm.** `NumberForm` is a **print wrapper**: `builtin_numberform`
(`src/numberform.c`) is inert — it returns `NULL` so the `NumberForm[...]` head
survives in the tree (`FullForm[NumberForm[1.23, 2]]` is itself), and all the work
happens at print time. `print.c` installs an active `NumberFormCtx` built from the
precision spec and option rules (`build_ctx` / `nf_parse_spec` /
`nf_apply_option_rule`) and routes every numeric leaf through
`numberform_render_number`. Per number, `nf_format_parts`: rejects non-finite,
special-cases zero and exact integers; extracts `count` significant base-10 digits
plus a decimal exponent (`mpfr_get_str` for MPFR, `"%.*e"` for machine reals);
chooses scientific vs decimal (`ScientificNotationThreshold` / `ExponentFunction`
/ `ExponentStep`); lays the digits into integer and fractional strings (a `{n, f}`
spec re-rounds to `f` fractional digits, a plain `n` drops trailing zeros); then
applies `DigitBlock` grouping and assembles sign / `NumberPoint` /
`NumberMultiplier` / `NumberFormat`.

**Data structures.** A small growable string buffer (`SB`) builds each rendered
number; the `NumberFormCtx` holds the parsed spec, the option strings, and the
measured alignment field. The *same* `nf_format_parts` drives both the measure
pass (`nf_measure`, for padding/alignment widths) and the output pass, so the
printed widths can never disagree.

**Complexity / limits.** Linear in the digit count per number and in the leaf
count of the wrapped expression. A requested precision below the integer-digit
count issues `NumberForm::reqsigz` (once per print, via `mth_message`) and pads
with zeros. Because the head survives, an intervening `NumberForm` blocks
arithmetic on the surrounding expression. Attributes `NHoldRest`, `Protected` —
the first argument evaluates, the spec and options are held under numeric
evaluation.

- `NHoldRest`, `Protected`. The first argument evaluates; the precision spec and options are held under numeric evaluation.
- A requested precision lower than the integer-digit count issues `NumberForm::reqsigz` and pads with zeros (`NumberForm[12345.6, 3]` is `12300.`).

**Attributes:** `NHoldRest`, `Protected`.

## References

- Source: [`src/numberform.c`](https://github.com/stblake/mathilda/blob/main/src/numberform.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_numberform.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numberform.c)
- Tests: [`tests/test_print.c`](https://github.com/stblake/mathilda/blob/main/tests/test_print.c)

## Notes & additional examples

### Notes

`NumberForm[expr, n]` displays the approximate reals in `expr` to `n` significant
digits; `NumberForm[expr, {n, f}]` uses `n` significant figures shown with exactly
`f` digits after the point. It works over integers, scalars, lists, matrices, and
mixed symbolic expressions — every inexact real inside `expr` is reformatted — and
takes a rich option set (`DigitBlock`, `NumberSeparator`, `NumberFormat`,
`ScientificNotationThreshold`, and more).

It is an inert **print wrapper**: the `NumberForm[...]` head stays in the
expression tree (so `FullForm` shows it) and only changes how the wrapped value is
displayed. Because the head survives, an intervening `NumberForm` blocks arithmetic
on the surrounding expression, so assign a variable first if the result must stay
computable. A requested precision below the integer-digit count issues
`NumberForm::reqsigz` and pads with zeros.
