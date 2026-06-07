### Worked examples

```mathematica
In[1]:= TeXForm[x^2 + 1]
Out[1]= 1+x^{2}

In[2]:= TeXForm[1/2]
Out[2]= \frac{1}{2}

In[3]:= TeXForm[Sqrt[x]/y]
Out[3]= \frac{\sqrt{x}}{y}
```

### Notes

`TeXForm` renders an expression as a TeX/LaTeX string, mapping fractions to `\frac`, roots to `\sqrt`, and superscripts to `^{...}`. Paste the output directly into a math environment.
