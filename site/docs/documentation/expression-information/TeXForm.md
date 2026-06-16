# TeXForm

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
TeXForm[expr]
    prints as a TeX version of expr.
TeXForm produces AMS-LaTeX-compatible TeX output.
When an input evaluates to TeXForm[expr], TeXForm does not appear in the output.
TeXForm translates standard mathematical functions and operations.
Following standard mathematical conventions, single-character symbol names are given in italic font, while multiple character names are given in roman font.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `TeXForm[expr]` is a held wrapper: `builtin_texform` returns `NULL`, and the actual rendering is triggered by `print_standard` in `src/print.c`, which on seeing the `TeXForm` head dispatches to the recursive `print_tex` renderer. `print_tex` emits AMS-LaTeX, reusing the same precedence levels as `get_expr_prec` so parentheses are inserted consistently with the ordinary printer. Numbers go through `print_tex_number` (integers/bigints verbatim, reals via `%g` with a trailing dot, MPFR at precision-derived digit counts). Symbols go through `print_tex_symbol`, which maps named constants (`Pi -> \pi`, `Infinity -> \infty`, `EulerGamma -> \gamma`, `Degree -> {}^{\circ}`, booleans to `\text{...}`, etc.), emits single-character names bare (italicised by math mode) and wraps multi-character names in `\text{...}`. Trig and hyperbolic heads are mapped to `\sin`, `\cosh`, … (with inverse variants rendered as `\sin^{-1}` style) via the `tex_trig_lookup` table; structural forms (`Power`, `Times`, fractions, `Sqrt`, lists/matrices, `SeriesData`) have dedicated emitters. The output is printed directly to stdout rather than returned as a string expression.

**Attributes:** `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/print.c`](https://github.com/stblake/mathilda/blob/main/src/print.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= TeXForm[x^2 + 1]
Out[1]= 1+x^{2}

In[2]:= TeXForm[1/2]
Out[2]= \frac{1}{2}

In[3]:= TeXForm[Sqrt[x]/y]
Out[3]= \frac{\sqrt{x}}{y}
```

```mathematica
In[1]:= TeXForm[Integrate[1/(1 + x^2), x]]
Out[1]= \tan ^{-1}(x)
```

```mathematica
In[1]:= TeXForm[Gamma[1/2]]
Out[1]= \sqrt{\pi }
```

```mathematica
In[1]:= TeXForm[D[Tan[x], x]]
Out[1]= \sec(x)^{2}
```

### Notes

`TeXForm` renders an expression as a TeX/LaTeX string, mapping fractions to `\frac`, roots to `\sqrt`, and superscripts to `^{...}`. Paste the output directly into a math environment.
