(* Experiment 25 -- Vector math for the transcendentals.
   ROADMAP ITEM 8 -- "a vector math library for the transcendentals: Accelerate's
   vv*, already linked", valued at ~1.9x on every elementwise Sin/Exp/Log in the
   system (experiments 2 and 6).

   The claim is narrow and checkable: Accelerate is ALREADY linked (the build
   passes -framework Accelerate), and it ships vvsin/vvexp/vvlog which are
   roughly twice the throughput of a scalar libm loop.  If these rows sit at ~2x
   behind NumPy -- which does call a vector library -- the item is confirmed and
   its size is measured. *)

Get["../harness.m"];
Get["../data.m"];

require[{"Sin", "Cos", "Tan", "Exp", "Log", "Sqrt", "ArcTan", "Sinh", "Power"}];

n = 4000000;
v = rand01[{n}] + 1/2;

bench["Sin over 4x10^6", Sin[v];];
check["Sin over 4x10^6", Round[10^6 N[Sin[1/2]]]];

bench["Exp over 4x10^6", Exp[v];];
check["Exp over 4x10^6", Round[10^6 N[Exp[1/2]]]];

bench["Log over 4x10^6", Log[v];];
check["Log over 4x10^6", Round[10^6 N[Log[3/2]]]];

bench["Sqrt over 4x10^6", Sqrt[v];];
check["Sqrt over 4x10^6", Round[10^6 N[Sqrt[3/2]]]];

bench["Tan over 4x10^6", Tan[v];];
check["Tan over 4x10^6", Round[10^6 N[Tan[1/2]]]];

bench["ArcTan 2-arg over 4x10^6", ArcTan[v, v + 1];];
check["ArcTan 2-arg over 4x10^6", Round[10^6 N[ArcTan[1/2, 3/2]]]];

(* A composite: three transcendentals in one expression.  If the system fuses,
   this costs one pass; if not, three passes plus two temporaries. *)
bench["Sin[Exp[Log[v]]] fused?", Sin[Exp[Log[v]]];];
check["Sin[Exp[Log[v]]] fused?", Round[10^6 N[Sin[Exp[Log[3/2]]]]]];

(* Power with a non-integer exponent: the pow() path, not repeated multiply. *)
bench["v^2.5 over 4x10^6", v^2.5;];
check["v^2.5 over 4x10^6", Round[10^6 N[(3/2)^2.5]]];
