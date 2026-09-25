# Higher-level mathematics ships as rewrite rules in Mathilda's own language.
# The derivative rules in src/internal/deriv.m are ordinary DownValues:
D[Sin[x], x]
D[x^n, x]
D[f[g[x]], x]
# The integral tables likewise load from .m files at startup:
Integrate[1/(1 + x^2), x]
