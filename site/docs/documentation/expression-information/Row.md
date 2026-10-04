# Row

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Row[{e1, e2, ...}] displays the ei concatenated together in a row.`**

Row\[{e1, e2, ...}, s\] inserts the string s between successive elements.  Strings are shown without quotes.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Row[{"x", " = ", 5}]
Out[1]= x = 5

In[2]:= Row[{1, 2, 3}, ", "]
Out[2]= 1, 2, 3
```

### Applications (3)

Elements concatenated with no separator

```mathematica
In[3]:= Row[{a, b, c}]
Out[3]= abc
```

A string inserted between successive elements

```mathematica
In[4]:= Row[{1, 2, 3}, ", "]
Out[4]= 1, 2, 3
```

Strings print without their quotes

```mathematica
In[5]:= Row[{"x", "=", 5}]
Out[5]= x=5
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

**Algorithm.** `Row` is a display wrapper like `NumberForm`: `builtin_row`
(`src/numberform.c`) is inert (returns `NULL`, so the head stays in the tree) and
the rendering lives in the standard printer. When `print.c` meets
`Row[{e1, e2, ...}]` (first argument a `List`, and only outside `InputForm`) it
sets the output-form flag and prints each element with `print_standard`, so
strings appear **without quotes**; a second string argument `Row[{...}, s]` is
emitted between successive elements as a separator. Outside that shape the literal
`Row[...]` prints normally.

**Data structures.** No allocation of its own — it writes directly to `stdout`
through the printer, toggling and restoring the global `g_print_output_form`
flag around the element loop. Its chief internal use is as the assembler for
`NumberForm`'s `NumberFormat` option, where `numberform_format_result_to_string`
renders a `Row` of mantissa / base / exponent pieces to an OutputForm string.

**Complexity / limits.** Linear in the number of elements. `Row` is a pure
presentation head: it has no computational value and takes the list's first
argument as its content, with an optional string separator. Attributes
`Protected`.

**Attributes:** `Protected`.

## References

**See also:** [NumberForm](../../expression-information/NumberForm/)

- Source: [`src/numberform.c`](https://github.com/stblake/mathilda/blob/main/src/numberform.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_numberform.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numberform.c)

## Notes & additional examples

### Notes

`Row[{e1, e2, ...}]` displays the elements concatenated left to right, with
strings shown without quotes; `Row[{...}, s]` inserts the string `s` between
successive elements. It is a pure presentation head with no computational value,
and its main internal use is assembling a custom display for `NumberForm`'s
`NumberFormat` option — the pieces `(mantissa, "10", exponent)` are handed to the
`NumberFormat` function, which typically wraps them in a `Row`.
