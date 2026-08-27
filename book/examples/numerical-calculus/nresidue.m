# Section 4.5.7: NResidue -- residues by a small contour integral.
Chop[NResidue[1/x, {x, 0}]]
Chop[NResidue[Tan[z], {z, Pi/2}]]
Chop[NResidue[Exp[1/x], {x, 0}, Radius -> 1]]
Chop[NResidue[Exp[x]/x^4, {x, 0}, Radius -> 1]]
