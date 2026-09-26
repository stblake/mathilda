# 4.2 Refine: simplifying under assumptions
Refine[Sqrt[x^2], x > 0]
Refine[Sqrt[x^2], Element[x, Reals]]
Refine[Abs[x], x < 0]
Refine[Log[x y], x > 0]
Assuming[x > 0, Refine[Sqrt[x^2 y^2], y < 0]]
