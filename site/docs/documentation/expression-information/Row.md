# Row

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Row[{e1, e2, ...}] displays the ei concatenated together in a row.`**

Row\[{e1, e2, ...}, s\] inserts the string s between successive elements.  Strings are shown without quotes.

## Examples

_No verified examples yet for this function._

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

**Attributes:** `Protected`.

## References

**See also:** [NumberForm](../../expression-information/NumberForm/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_numberform.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numberform.c)
