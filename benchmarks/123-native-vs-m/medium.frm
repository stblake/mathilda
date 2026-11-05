* medium.frm -- Rolf Mertig, Woxi issue #2142 (reference only; run_all.py
* has no FORM runner).   time form -q medium.frm
* Naive recursive Dirac trace, same algorithm as the Wolfram-language diracTr.
* Tr() = 4
* Tr(a) = 0
* Tr(a, b1..bn) = sum_m (-1)^(m+1) * g(a,bm) * Tr(without bm)

Symbols x1,...,x14, a,b, END;
CFunctions g, Tr, Wp, Wm;

Local res = Tr(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14);

repeat;
  id  Tr = 4;
  id  Tr(a?,?r) = Wp(a,?r,END);

* current sign +, then continue with -
  id  Wp(a?,b?,?r,END,?d) =  g(a,b)*Tr(?d,?r) + Wm(a,?r,END,?d,b);
  id  Wp(a?,END,?d) = 0;

* current sign -, then continue with +
  id  Wm(a?,b?,?r,END,?d) = -g(a,b)*Tr(?d,?r) + Wp(a,?r,END,?d,b);
  id  Wm(a?,END,?d) = 0;
endrepeat;

.sort
.end
