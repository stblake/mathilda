# NumberForm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NumberForm[expr, n] prints approximate real numbers in expr to n-digit`**

precision; NumberForm\[expr, {n, f}\] uses n digits with f to the right of the decimal point; NumberForm\[expr\] uses the default options. Works on integers as well.  It is an inert print wrapper: the head remains in the expression and only changes how it is displayed. Options: DefaultPrintPrecision, DigitBlock, ExponentFunction, ExponentStep, NumberFormat, NumberMultiplier, NumberPadding, NumberPoint, NumberSeparator, NumberSigns, ScientificNotationThreshold, SignPadding.

## Examples (3)

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

- `NHoldRest`, `Protected`. The first argument evaluates; the precision spec and options are held under numeric evaluation.
- A requested precision lower than the integer-digit count issues `NumberForm::reqsigz` and pads with zeros (`NumberForm[12345.6, 3]` is `12300.`).

**Attributes:** `NHoldRest`, `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_numberform.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numberform.c)
