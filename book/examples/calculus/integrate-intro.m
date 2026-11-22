# 5.3.5 Indefinite integration -- a first look; the full story is the Symbolic Integration chapter
Integrate[1/(x^4 + x^2 + 1), x]
Integrate[Exp[-x^2], x]
Integrate[1/(x^4 + 1)^(1/4), x]
Check[Integrate[Sqrt[x + Log[x]], x], Last[$MessageList]]
