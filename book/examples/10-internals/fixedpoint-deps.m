# Dependency-recorded fixed points: a cached verdict is dropped exactly when what it read changes
x = 5;
Block[{x}, x + y]
k[z_] := 0 /; !TrueQ[flag]
flag = True;
k[2]
Block[{flag}, k[2]]
k[2]
flag = False;
k[2]
