# InputForm prints machine reals with the shortest digits that read back exactly.
0.1 + 0.2
InputForm[0.1 + 0.2]
InputForm[{0.1, 1/3., 1.5*^-5}]
ToExpression[ToString[0.1 + 0.2, InputForm]] === 0.1 + 0.2
# Up to 17 significant digits a literal is a machine real; 18 or more make it arbitrary precision (v0.6136):
Precision[0.30000000000000004]
Precision[0.300000000000000044]
# A machine real's display re-reads as a Real (v0.6604):
1.*^15
