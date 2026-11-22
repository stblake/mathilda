# Ch.7 -- packages and contexts: write a package to a temporary file, then load it
$Context
$ContextPath
Context[Sin]
Context[NDArray]
Context[SeparablePolynomialDecomposition]
MemberQ[$ContextPath, "SeparablePolynomialDecomposition`"]
pkg = FileNameJoin[{$TemporaryDirectory, "Geometry.m"}];
WriteString[pkg, "BeginPackage[\"Geometry`\"]\nDiskArea::usage = \"DiskArea[r] is the area of a disk of radius r.\";\nBegin[\"`Private`\"]\nsquare[x_] := x^2\nDiskArea[r_] := Pi square[r]\nEnd[]\nEndPackage[]\n"]; Close[pkg];
Get[pkg]
DiskArea[3]
DiskArea::usage
Context[DiskArea]
First[$ContextPath]
MemberQ[$Packages, "Geometry`"]
square[3]
Geometry`Private`square[3]
Names["Geometry`*"]
Needs["Geometry`"]
$Context
FullForm[Hold[<< Geometry.m]]
Needs["Lean`"]
Take[$ContextPath, 2]
