# The printer chooses parentheses from precedence; the forms differ in detail.
InputForm[a/b]
FullForm[a/b]
(a + b) c
InputForm[(a + b) c]
FullForm[Rational[3, 4]]
3/4
FullForm[Complex[2, 3]]
2 + 3 I
InputForm[{1, 2, 3}]
