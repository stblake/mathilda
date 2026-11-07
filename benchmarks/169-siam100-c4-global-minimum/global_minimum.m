(* bench-timeout: 400 *)
(* ==========================================================================
   Experiment 169 -- SIAM 100-Digit Challenge, Problem 4 ("Think Globally,
   Act Locally").

       min over (x, y) of  e^sin(50x) + sin(60 e^y) + sin(70 sin x)
                           + sin(sin(80y)) - sin(10(x+y)) + (x^2+y^2)/4
                                                       =  -3.3068686474...

   Idiomatic statement: the system's GLOBAL minimiser on the box [-1,1]^2
   (WL NMinimize with the box as constraints; scipy differential_evolution
   with bounds).  The box is the user's one piece of insight: outside it
   (x^2+y^2)/4 > 1/4 while the other terms are >= -4 + 1/e, and the global
   minimum lies well inside (at (-0.0244, 0.2106)); the objective has
   ~2700 local minima in the box.
     c4.machine   NMinimize at machine precision;
     c4.wp50      NMinimize with WorkingPrecision -> 60.
   Maxima and FriCAS have no global minimiser: SKIP.  mpmath has no minimiser
   either; its wp50 row is the machine global search followed by
   mpmath.findroot on the gradient (mpmath.diff), the WorkingPrecision
   polish WL applies internally.
   True value and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NMinimize"}];

f4 = Exp[Sin[50 x]] + Sin[60 Exp[y]] + Sin[70 Sin[x]] + Sin[Sin[80 y]] -
     Sin[10 (x + y)] + (x^2 + y^2)/4;

vrun["c4.machine", First[NMinimize[{f4, -1 <= x <= 1, -1 <= y <= 1}, {x, y}]], 10, 3];
vrun["c4.wp50", First[NMinimize[{f4, -1 <= x <= 1, -1 <= y <= 1}, {x, y},
                                WorkingPrecision -> 60]], 48];
