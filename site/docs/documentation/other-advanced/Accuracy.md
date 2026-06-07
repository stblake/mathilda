# Accuracy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Accuracy[x]
    Returns the number of digits of accuracy in x — equal to
    Precision[x] − Log10[Abs[x]]. Exact numbers (including exact 0)
    return Infinity. Inexact zeros are finite: machine 0. returns
    ≈ 323.607; MPFR 0 of precision p returns p digits.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_accuracy` (`src/precision.c`) delegates to `accuracy_of`, which returns the number of correct digits to the right of the decimal point. Exact quantities (integers, bigints, exact rationals, strings, symbols, exact zero) return `Infinity`. For inexact reals it returns `MachinePrecisionDigits - log10|x|`; for `EXPR_MPFR` it uses the value's actual precision (`mpfr_get_prec / log2(10)`) minus `log10|x|` (inexact zero gets the precision in digits directly). `Complex[re, im]` and general function arguments recurse and take the minimum (`precision_min`) over components/arguments. `Accuracy` is `ATTR_LISTABLE`, so it threads over lists.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/precision.c`](https://github.com/stblake/mathilda/blob/main/src/precision.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Accuracy[1.5]
Out[1]= 15.7785

In[2]:= Accuracy[5]
Out[2]= Infinity

In[3]:= Accuracy[N[Pi, 30]]
Out[3]= 29.6058
```

### Notes

`Accuracy` gives the number of digits to the right of the decimal point, i.e. the digit count relative to the absolute (not relative) magnitude. Exact numbers have `Infinity` accuracy; for a fixed precision, accuracy decreases as the magnitude grows.
