# The numeric tower after the 2026-11 review.
1.5`20
InputForm[1.5`20]
1.*^15
x = 1.5`20 - 1.5`20
Accuracy[x]
y = Exp[1000.]
Precision[y]
MachineNumberQ[y]
MachineNumberQ[Exp[1000.] Exp[-999.]]
N[Cos[10^-20] - 1]
N[Cos[10^-20] - 1, 20]
Exp[10.^10]
Exp[-10.^10]
Precision[{1.5`20, 2.5`30}]
Block[{$MaxPrecision = 30}, Precision[N[Pi, 50]]]
