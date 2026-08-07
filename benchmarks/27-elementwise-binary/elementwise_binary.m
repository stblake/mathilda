(* Experiment 27 -- Elementwise binary kernels.
   ROADMAP ITEM 10 -- "vectorise the elementwise binary kernels: MapThread[Min]
   runs at 7.7 GB/s where NumPy reaches 48", valued at 6.4x (experiment 19).

   Binary elementwise is where a system's SIMD story shows up most plainly: the
   loop body is one instruction, so anything other than a vector op is visible
   immediately.  Min/Max in particular have no arithmetic at all -- they are a
   compare and a select, which vectorises perfectly. *)

Get["../harness.m"];
Get["../data.m"];

require[{"MapThread", "Min", "Max", "Plus", "Times", "Mod", "Quotient",
         "Clip", "UnitStep", "Threshold"}];

n = 4000000;
a = rand01[{n}];
b = rand01[{n}];

bench["a + b over 4x10^6", a + b;];
check["a + b over 4x10^6", Total[Range[5] + Range[5]]];

bench["a * b over 4x10^6", a b;];
check["a * b over 4x10^6", Total[Range[5] Range[5]]];

bench["MapThread[Min] over 4x10^6", MapThread[Min, {a, b}];];
check["MapThread[Min] over 4x10^6",
  Total[MapThread[Min, {{1, 5, 3}, {4, 2, 6}}]]];

bench["MapThread[Max] over 4x10^6", MapThread[Max, {a, b}];];
check["MapThread[Max] over 4x10^6",
  Total[MapThread[Max, {{1, 5, 3}, {4, 2, 6}}]]];

(* Clip is the same compare-and-select against constants rather than an array.
   The bounds are MACHINE REALS, deliberately. Written as the exact {1/4, 3/4}
   this row measured something else entirely: an exact bound is preserved at
   every clipped position, so Clip[reals, {1/4,3/4}] returns a MIXED list --
   {1/4, 0.5, 3/4}, Rationals intact -- which by definition cannot be a packed
   buffer, while the numpy column clips float64 to float64. That is the same
   exact-vs-float64 mismatch experiment 20 documents for Fit, and it cost 700x:
   0.193 s against 0.000256 s for the identical clamp with {0.25, 0.75}, where
   Mathilda is ~3.7x FASTER than numpy. The exact-bounds cost is real but it is
   the price of exactness, not a kernel defect, and it is not what this
   experiment measures. *)
bench["Clip to [0.25, 0.75] over 4x10^6", Clip[a, {0.25, 0.75}];];
check["Clip to [0.25, 0.75] over 4x10^6",
  Round[10^6 N[Total[Clip[{0., 1/2, 1.}, {0.25, 0.75}]]]]];

(* A three-operand expression: tests whether temporaries are avoided. *)
bench["a b + a over 4x10^6", a b + a;];
check["a b + a over 4x10^6", Total[Range[5] Range[5] + Range[5]]];

(* Integer elementwise: the int64 buffer rather than float64. *)
ia = RandomInteger[{1, 1000}, {n}];
ib = RandomInteger[{1, 1000}, {n}];
bench["integer Mod over 4x10^6", Mod[ia, ib];];
check["integer Mod over 4x10^6", Total[Mod[{7, 8, 9}, {3, 5, 4}]]];
