(* ==========================================================================
   Experiment 56 -- Multidimensional & hard quadrature (group E)
   ==========================================================================
   WHAT IT MEASURES.  NIntegrate beyond the 1-D smooth case of experiment 12:
   the 2-D/3-D cubature path (src/numerical_calculus/cubature.c), dependent
   integration bounds (a triangle and a disk), an oscillatory integrand
   (oscint.c / levincoll.c), an integrable endpoint singularity (dequad.c
   tanh-sinh), and a semi-infinite range (denint.c).

   Baseline is scipy.integrate (compiled quadrature), but its multidimensional
   drivers dblquad/tplquad call a PYTHON integrand once per node, so on the
   cubature rows this is partly a compiled-integrand vs interpreted-integrand
   comparison -- read the ratio with that in mind; the 1-D rows are the cleaner
   compiled-vs-compiled ones.

   Checks are Round[10^6 value]; every integrand has a known closed value, so
   the three systems pin the same number.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"NIntegrate"}];

bench["NIntegrate 2D x y on unit square", NIntegrate[x y, {x, 0, 1}, {y, 0, 1}];];
check["NIntegrate 2D x y on unit square", Round[10^6 NIntegrate[x y, {x, 0, 1}, {y, 0, 1}]]];

bench["NIntegrate 3D x y z on unit cube",
      NIntegrate[x y z, {x, 0, 1}, {y, 0, 1}, {z, 0, 1}];];
check["NIntegrate 3D x y z on unit cube",
      Round[10^6 NIntegrate[x y z, {x, 0, 1}, {y, 0, 1}, {z, 0, 1}]]];

bench["NIntegrate 2D Gaussian on [-3,3]^2",
      NIntegrate[Exp[-(x^2 + y^2)], {x, -3, 3}, {y, -3, 3}];];
check["NIntegrate 2D Gaussian on [-3,3]^2",
      Round[10^6 NIntegrate[Exp[-(x^2 + y^2)], {x, -3, 3}, {y, -3, 3}]]];

bench["NIntegrate triangle (dependent bounds)",
      NIntegrate[1, {x, 0, 1}, {y, 0, x}];];
check["NIntegrate triangle (dependent bounds)",
      Round[10^6 NIntegrate[1, {x, 0, 1}, {y, 0, x}]]];

bench["NIntegrate unit disk (dependent bounds)",
      NIntegrate[1, {x, -1, 1}, {y, -Sqrt[1 - x^2], Sqrt[1 - x^2]}];];
check["NIntegrate unit disk (dependent bounds)",
      Round[10^6 NIntegrate[1, {x, -1, 1}, {y, -Sqrt[1 - x^2], Sqrt[1 - x^2]}]]];

bench["NIntegrate oscillatory Sin[50x]/x on [1,10]",
      NIntegrate[Sin[50 x]/x, {x, 1, 10}];];
check["NIntegrate oscillatory Sin[50x]/x on [1,10]",
      Round[10^6 NIntegrate[Sin[50 x]/x, {x, 1, 10}]]];

bench["NIntegrate singular 1/Sqrt[x] on [0,1]",
      NIntegrate[1/Sqrt[x], {x, 0, 1}];];
check["NIntegrate singular 1/Sqrt[x] on [0,1]",
      Round[10^6 NIntegrate[1/Sqrt[x], {x, 0, 1}]]];

bench["NIntegrate semi-infinite Exp[-x^2] on [0,Inf]",
      NIntegrate[Exp[-x^2], {x, 0, Infinity}];];
check["NIntegrate semi-infinite Exp[-x^2] on [0,Inf]",
      Round[10^6 NIntegrate[Exp[-x^2], {x, 0, Infinity}]]];
