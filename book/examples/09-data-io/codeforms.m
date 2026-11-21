# FortranForm and CForm: an expression as source code for another language.
FortranForm[x^2 + Sin[y]/3]
CForm[x^2 + Sin[y]/3]
FortranForm[Sqrt[1 - x^2] + Exp[-x]]
CForm[Mod[n, 7] + v[[3]] && a == b]
FortranForm[(1 + 2 I) z]
# ToString with FormatType gives the code as a string, ready to write to a file.
ToString[D[ArcTan[x/y], x], FormatType -> CForm]
