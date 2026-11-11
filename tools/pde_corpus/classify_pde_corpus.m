(* classify_pde_corpus.m -- structural class of every PDE corpus record.
 *
 *   ./Mathilda -file tools/pde_corpus/classify_pde_corpus.m > DSolve_test_status/pde/classes.tsv
 *
 * One TSV line per record: label, class, nvars, order, ncond, domain.  No DSolve
 * call: the class is read off the equation, so it is cheap and deterministic,
 * and it buckets the corpus by the METHOD a record needs (DSOLVE_PDE_PLAN.md
 * waves), which tools/dsolve_pde_corpus_report.py joins with the verdicts.
 *
 * class, for a single PDE:
 *   o1-lin-const / o1-lin-var    first order, linear, constant / variable coeffs
 *   o1-quasi                     first order, linear in u_v's, coefficients with u
 *   o1-nonlin                    first order, nonlinear in the derivatives
 *   o2-lin-const-hom / -inh      second order, linear, constant coefficients
 *   o2-lin-var                   second order, linear, variable coefficients
 *   o2-nonlin                    second order, nonlinear
 *   oN-lin / oN-nonlin           order >= 3
 *   system                       more than one PDE
 * domain: none (general solution), ivp (data on one line), bounded (two lines of
 * one variable: an interval, rectangle, ... -> eigenfunction expansions), other.
 *)
pcFiles = FileNames["DSolve_test_status/pde/DE_PDE_*.m"];

pcU[fn_Symbol] := fn;  pcU[fn_] := Head[fn];
pcIsCond[eq_, u_, vars_] := FreeQ[eq, u @@ vars] && FreeQ[eq, Derivative[__][u] @@ vars];

pcClass[eqns_, fn_, vars_] := Module[
  {eqs = Flatten[{eqns}], u = pcU[fn], pdes, conds, terms, syms, R, ord, hess, lin,
   d1, dsym, quasi, const, hom, slots, dom, cls},
  pdes  = Select[eqs, !pcIsCond[#, u, vars] &];
  conds = Select[eqs, pcIsCond[#, u, vars] &];
  ord = Max[0, Cases[pdes, Derivative[n__][u][___] :> Total[{n}], {0, Infinity}]];
  (* domain: which slots the conditions fix, and at how many values *)
  slots = Cases[conds, (u | Derivative[__][u])[a__] /; Length[{a}] == Length[vars] :>
            Select[Range[Length[vars]], {a}[[#]] =!= vars[[#]] &] -> {a}, {0, Infinity}];
  dom = Which[
    conds === {}, "none",
    AnyTrue[Range[Length[vars]], Function[j,
      Length[Union@Cases[slots, ({j} -> a_) :> a[[j]]]] >= 2]], "bounded",
    Length[Union[First /@ slots]] == 1, "ivp",
    True, "other"];
  If[Length[pdes] != 1,
    Return[{If[Length[pdes] > 1, "system", "none"], Length[vars], ord, Length[conds], dom}]];
  R = pdes[[1]] /. Equal -> Subtract;
  terms = Union@Cases[R, (u | Derivative[__][u])[Sequence @@ vars], {0, Infinity}];
  syms = Table[Unique["pc"], {Length[terms]}];
  R = R /. Thread[terms -> syms];
  If[!FreeQ[R, u], Return[{"nonlocal", Length[vars], ord, Length[conds], dom}]];
  hess = Quiet@Flatten@Table[Together[D[R, a, b]], {a, syms}, {b, syms}];
  lin = AllTrue[hess, # === 0 &];
  const = lin && AllTrue[D[R, #] & /@ syms, FreeQ[#, Alternatives @@ vars] &];
  hom = lin && Together[R /. Thread[syms -> 0]] === 0;
  (* quasilinear: affine in the first-derivative symbols (u itself may appear) *)
  dsym = Pick[syms, (FreeQ[#, u @@ vars] &) /@ terms];
  quasi = AllTrue[Quiet@Flatten@Table[Together[D[R, a, b]], {a, dsym}, {b, dsym}], # === 0 &];
  cls = Which[
    ord == 1 && lin && const, "o1-lin-const",
    ord == 1 && lin,          "o1-lin-var",
    ord == 1 && quasi,        "o1-quasi",
    ord == 1,                 "o1-nonlin",
    ord == 2 && const && hom, "o2-lin-const-hom",
    ord == 2 && const,        "o2-lin-const-inh",
    ord == 2 && lin,          "o2-lin-var",
    ord == 2,                 "o2-nonlin",
    lin,                      "oN-lin",
    True,                     "oN-nonlin"];
  {cls, Length[vars], ord, Length[conds], dom}
];

(* Iterators are pc-prefixed: a corpus record spells r (Black-Scholes's rate),
 * f, c, ... as ITS symbols, so `Do[..., {r, recs}]` would bind r to a record
 * containing r and recurse without end. *)
Do[Do[With[{pcC = Quiet@TimeConstrained[pcClass[pcRec[[3]], pcRec[[4]], pcRec[[5]]], 10,
                                         {"timeout", Length[pcRec[[5]]], "?", "?", "?"}]},
     Print[StringRiffle[ToString /@ Prepend[pcC, pcRec[[1]]], "\t"]]],
   {pcRec, Get[pcFile]}], {pcFile, pcFiles}];
