# Attributes are the small vocabulary that steers the generic evaluator.
Attributes[Plus]
Attributes[Sin]
Attributes[Table]
Attributes[Set]
# Orderless sorts arguments into canonical order:
b + a + c
# Flat flattens nested same-head calls -- before, then after opting in:
f[f[a, b], c]
SetAttributes[f, Flat]
f[f[a, b], c]
# Listable threads a head over list arguments:
h[{1, 2, 3}]
SetAttributes[h, Listable]
h[{1, 2, 3}]
