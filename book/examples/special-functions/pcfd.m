# Parabolic cylinder functions D_nu(z)
ParabolicCylinderD[nu, 0]
ParabolicCylinderD[2, z]
Simplify[FunctionExpand[ParabolicCylinderD[2, z]]]
FunctionExpand[ParabolicCylinderD[-1, z]]
ParabolicCylinderD[2.5, 3.]
N[ParabolicCylinderD[1/3, 2], 30]
D[ParabolicCylinderD[nu, z], z]
# Outside the double range the value is still returned
ParabolicCylinderD[2.5, 60.]
