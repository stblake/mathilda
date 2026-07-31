#!/usr/bin/env python3
"""
Mathilda vs Mathematica vs NumPy on high-performance-computing kernels.

Runs each benchmark in all three systems, one process per system, timing with
each system's wall clock -- AbsoluteTiming in the two CAS, time.perf_counter in
Python. NOT Timing[], which reports CPU time summed over threads and so
over-reports every threaded NDArray path and every BLAS call by roughly the core
count.

Mathilda and Mathematica run the SAME source text for every benchmark; a `wl`
override is used only where the two languages genuinely differ. The NumPy column
is a separate implementation of the same algorithm -- same operations in the same
order, so the comparison is of execution, not of method. Where a benchmark has a
cheap scalar answer that answer is recorded from every system and compared; a
timing row is only meaningful if they agree.

Why three systems: Mathematica says whether we are behind a competitor. NumPy --
which on this host links the same Accelerate BLAS Mathilda does -- says whether
we are behind the machine. On the dense linear algebra rows all three call
identical kernels, so any spread there is pure overhead; on the array rows NumPy
is a thin memory-bandwidth reference.

Usage:
    comparisons/hpc_bench.py                 # full run
    comparisons/hpc_bench.py --scale 0.1     # smaller sizes, for a smoke test
    comparisons/hpc_bench.py --only fft,sort # a subset, by id
    comparisons/hpc_bench.py --system mathilda,numpy
Output: a markdown table on stdout, and the raw JSON to --json.
"""
import argparse, json, os, re, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MATHILDA = os.environ.get("MATHILDA_BIN", os.path.join(ROOT, "Mathilda"))
WOLFRAM = os.environ.get(
    "WOLFRAMSCRIPT", "/Applications/Mathematica.app/Contents/MacOS/wolframscript")
# A python that actually has numpy -- NOT necessarily the one running this script.
PYTHON = os.environ.get("HPC_PYTHON", "/usr/local/bin/python3.11")

# --------------------------------------------------------------------------
# Benchmarks.
#
# id        short key for --only
# group     section heading in the report
# name      row label
# setup     run once, untimed (data construction, kernel definition)
# expr      the timed expression; ends in ';' where the value is a big array
# check     optional: a cheap scalar recorded from every system and compared
# reps      timed repetitions; the reported number is the MINIMUM
# py        NumPy translation of `expr` -- a Python EXPRESSION (helpers go in
#           py_setup, exactly as the CAS versions put them in `setup`).
#           Omitted where NumPy has no equivalent; the column then reads '—'.
# py_setup  NumPy translation of `setup`; statements, exec'd once
# py_check  NumPy translation of `check`; omit to skip the value comparison for
#           this system only (used where the libraries implement genuinely
#           different algorithms, e.g. two different cubic interpolants)
# --------------------------------------------------------------------------
B = []


def bench(**kw):
    kw.setdefault("reps", 5)
    kw.setdefault("check", None)
    kw.setdefault("wl", None)
    kw.setdefault("wl_setup", None)
    kw.setdefault("py", None)
    kw.setdefault("py_setup", None)
    kw.setdefault("py_check", None)
    kw.setdefault("note", "")
    kw.setdefault("cold", False)
    B.append(kw)


LA = "Dense linear algebra (BLAS / LAPACK)"
SP = "Spectral"
AR = "Array and memory primitives (n = 10^7)"
ST = "Stencils and PDE relaxation"
SC = "Scalar kernels via Compile[]"
IN = "Integer and combinatorial"
BN = "Arbitrary precision (GMP / MPFR)"

# ---- dense linear algebra -------------------------------------------------
bench(id="matmul", group=LA, name="Matrix multiply, 1000x1000",
      setup="n=1000; A=RandomReal[{0,1},{n,n}]; Bm=RandomReal[{0,1},{n,n}];",
      expr="A . Bm;", reps=5,
      py_setup="n=1000\nA=np.random.rand(n,n)\nBm=np.random.rand(n,n)",
      py="A @ Bm")
bench(id="solve", group=LA, name="LinearSolve, 1000x1000",
      setup="n=1000; A=RandomReal[{0,1},{n,n}]; bv=RandomReal[{0,1},n];",
      expr="LinearSolve[A,bv];", reps=5,
      py_setup="n=1000\nA=np.random.rand(n,n)\nbv=np.random.rand(n)",
      py="np.linalg.solve(A,bv)")
bench(id="inverse", group=LA, name="Inverse, 500x500",
      setup="n=500; A=RandomReal[{0,1},{n,n}];", expr="Inverse[A];", reps=5,
      py_setup="n=500\nA=np.random.rand(n,n)", py="np.linalg.inv(A)")
bench(id="det", group=LA, name="Det, 500x500",
      setup="n=500; A=RandomReal[{0,1},{n,n}];", expr="Det[A];", reps=5,
      py_setup="n=500\nA=np.random.rand(n,n)", py="np.linalg.det(A)")
bench(id="qr", group=LA, name="QRDecomposition, 500x500",
      setup="n=500; A=RandomReal[{0,1},{n,n}];", expr="QRDecomposition[A];", reps=3,
      py_setup="n=500\nA=np.random.rand(n,n)", py="np.linalg.qr(A)")
bench(id="eigen", group=LA, name="Eigenvalues, 300x300 symmetric",
      setup="n=300; A=RandomReal[{0,1},{n,n}]; S=(A+Transpose[A])/2.;",
      expr="Eigenvalues[S];", reps=3,
      py_setup="n=300\nA=np.random.rand(n,n)\nS=(A+A.T)/2.0",
      py="np.linalg.eigvalsh(S)")
bench(id="svd", group=LA, name="SingularValueDecomposition, 300x300",
      setup="n=300; A=RandomReal[{0,1},{n,n}];",
      expr="SingularValueDecomposition[A];", reps=3,
      py_setup="n=300\nA=np.random.rand(n,n)", py="np.linalg.svd(A)")

# ---- spectral -------------------------------------------------------------
bench(id="fft", group=SP, name="Fourier, 2^20 reals",
      setup="v=RandomReal[{0,1},2^20];", expr="Fourier[v];", reps=5,
      py_setup="v=np.random.rand(2**20)", py="np.fft.fft(v)")
bench(id="fft2", group=SP, name="Fourier, 2^20 complex",
      setup="vz=RandomReal[{0,1},2^20] + I RandomReal[{0,1},2^20];",
      expr="Fourier[vz];", reps=5,
      py_setup="vz=np.random.rand(2**20)+1j*np.random.rand(2**20)",
      py="np.fft.fft(vz)")

# ---- array / memory primitives -------------------------------------------
ARSET = "n=10^7; x=RandomReal[{0,1},n]; y=RandomReal[{0,1},n];"
ARSETP = "n=10**7\nx=np.random.rand(n)\ny=np.random.rand(n)"
bench(id="triad", group=AR, name="STREAM triad, a = b + 3 c",
      setup=ARSET, expr="x + 3. y;", reps=5,
      py_setup=ARSETP, py="x + 3.0*y")
bench(id="total", group=AR, name="Total (reduction)", setup=ARSET,
      expr="Total[x];", reps=5, py_setup=ARSETP, py="x.sum()")
bench(id="accum", group=AR, name="Accumulate (prefix scan)", setup=ARSET,
      expr="Accumulate[x];", reps=5, py_setup=ARSETP, py="np.cumsum(x)")
bench(id="sort", group=AR, name="Sort", setup=ARSET, expr="Sort[x];", reps=3,
      py_setup=ARSETP, py="np.sort(x)")
bench(id="sin", group=AR, name="Sin (elementwise)", setup=ARSET,
      expr="Sin[x];", reps=5, py_setup=ARSETP, py="np.sin(x)")
bench(id="exp", group=AR, name="Exp (elementwise)", setup=ARSET,
      expr="Exp[x];", reps=5, py_setup=ARSETP, py="np.exp(x)")
bench(id="dot", group=AR, name="Dot (inner product)", setup=ARSET,
      expr="x . y;", reps=5, py_setup=ARSETP, py="x @ y")
# The structural family. Every one of these was on the marshalling path (one Expr
# per element) until this round of work; they are the operations array-style
# numerical code is actually written out of.
#
# Several have a NumPy spelling that returns a VIEW rather than a new array
# (x[::-1], reshape). A view is O(1) and would flatter NumPy against two systems
# that both materialise, so those rows force a copy -- the comparison is of the
# same work, not of the same syntax.
bench(id="rotate", group=AR, name="RotateLeft", setup=ARSET,
      expr="RotateLeft[x,3];", reps=5, py_setup=ARSETP, py="np.roll(x,-3)")
bench(id="reverse", group=AR, name="Reverse", setup=ARSET,
      expr="Reverse[x];", reps=5, py_setup=ARSETP, py="x[::-1].copy()")
bench(id="joinv", group=AR, name="Join (two 10^7 vectors)", setup=ARSET,
      expr="Join[x,y];", reps=3, py_setup=ARSETP, py="np.concatenate([x,y])")
bench(id="part2", group=AR, name="Partition, length 2", setup=ARSET,
      expr="Partition[x,2];", reps=3, py_setup=ARSETP, py="x.reshape(-1,2).copy()")
bench(id="diffs", group=AR, name="Differences", setup=ARSET,
      expr="Differences[x];", reps=3, py_setup=ARSETP, py="np.diff(x)")
bench(id="riffle", group=AR, name="Riffle with a Real", setup=ARSET,
      expr="Riffle[x,0.];", reps=3,
      py_setup=ARSETP + "\ndef riffle(a, f):\n"
                        "    r = np.full(2*a.size-1, f)\n"
                        "    r[0::2] = a\n    return r",
      py="riffle(x,0.0)")
bench(id="padr", group=AR, name="PadRight, Real fill", setup=ARSET,
      expr="PadRight[x,n+1,0.];", reps=3,
      py_setup=ARSETP, py="np.concatenate([x,np.zeros(1)])")
bench(id="padrd", group=AR, name="PadRight, default (exact 0) fill", setup=ARSET,
      expr="PadRight[x,n+1];", reps=3,
      note="mixed exact/inexact result: unpacked in BOTH systems, by definition")

# ---- stencils -------------------------------------------------------------
bench(id="jacobi", group=ST, name="Jacobi 5-point relaxation, 512^2, 100 sweeps",
      setup="n=512; u0=RandomReal[{0,1},{n,n}]; "
            "jac[u_]:=(RotateLeft[u,{1,0}]+RotateRight[u,{1,0}]"
            "+RotateLeft[u,{0,1}]+RotateRight[u,{0,1}])/4.;",
      expr="Nest[jac,u0,100];", reps=3,
      py_setup="n=512\nu0=np.random.rand(n,n)\n"
               "def jac(u):\n"
               "    return (np.roll(u,-1,0)+np.roll(u,1,0)"
               "+np.roll(u,-1,1)+np.roll(u,1,1))/4.0",
      py="nest(jac,u0,100)")
bench(id="life", group=ST, name="Game of Life, 256^2, 100 generations",
      setup="n=256; g0=RandomInteger[1,{n,n}]; "
            "nb[g_]:=Sum[RotateLeft[g,{i,j}],{i,-1,1},{j,-1,1}]-g; "
            "life[g_]:=With[{k=nb[g]}, "
            "UnitStep[k-3]UnitStep[3-k] + UnitStep[k-2]UnitStep[2-k] g];",
      expr="Nest[life,g0,100];", reps=3,
      py_setup="n=256\ng0=np.random.randint(0,2,(n,n))\n"
               "def nb(g):\n"
               "    s = -g\n"
               "    for i in (-1,0,1):\n"
               "        for j in (-1,0,1):\n"
               "            s = s + np.roll(np.roll(g,-i,0),-j,1)\n"
               "    return s\n"
               "def life(g):\n"
               "    k = nb(g)\n"
               "    return (k==3).astype(np.int64) + (k==2).astype(np.int64)*g",
      py="nest(life,g0,100)")

# ---- scalar kernels via Compile[] ----------------------------------------
# NumPy has no answer here and numba is not installed on this host, so the Python
# column is CPython running the same scalar loop. Those rows are labelled: they
# measure the interpreter, not the library, and reading them as a NumPy result
# would be wrong in both directions.
bench(id="logistic", group=SC, name="Logistic map, 10^7 iterations",
      setup="lg=Compile[{{x0,_Real},{m,_Integer}}, "
            "Module[{x=x0,k=0}, While[k<m, x=3.9 x (1.-x); k=k+1]; x]];",
      expr="lg[0.5, 10^7];", check="lg[0.5, 10^7]", reps=3,
      py_setup="def lg(x0, m):\n"
               "    x = x0\n"
               "    for _ in range(m): x = 3.9*x*(1.0-x)\n    return x",
      py="lg(0.5, 10**7)", py_check="lg(0.5, 10**7)",
      note="Python column is CPython, not NumPy: no numba on this host")
bench(id="mandel", group=SC, name="Mandelbrot, 800x800, 100 iterations",
      setup="mandel=Compile[{{cx,_Real},{cy,_Real},{mx,_Integer}}, "
            "Module[{zx=0.,zy=0.,t=0.,k=0}, "
            "While[k<mx && zx zx + zy zy < 4., "
            "t=zx zx - zy zy + cx; zy=2. zx zy + cy; zx=t; k=k+1]; k]]; "
            "n=800; h=2.5/(n-1);",
      expr="Table[mandel[x,yy,100],{yy,-1.25,1.25,h},{x,-2.,0.5,h}];",
      check="Total[Table[mandel[x,yy,100],{yy,-1.25,1.25,h},{x,-2.,0.5,h}],2]",
      reps=3,
      py_setup="n=800\nh=2.5/(n-1)\n"
               "def mandelgrid(m):\n"
               "    cy = (-1.25 + h*np.arange(n))[:,None]\n"
               "    cx = (-2.0 + h*np.arange(n))[None,:]\n"
               "    zx = np.zeros((n,n)); zy = np.zeros((n,n))\n"
               "    it = np.zeros((n,n), dtype=np.int64)\n"
               "    live = np.ones((n,n), dtype=bool)\n"
               "    for _ in range(m):\n"
               "        live &= (zx*zx + zy*zy) < 4.0\n"
               "        if not live.any(): break\n"
               "        t = zx*zx - zy*zy + cx\n"
               "        zy = np.where(live, 2.0*zx*zy + cy, zy)\n"
               "        zx = np.where(live, t, zx)\n"
               "        it += live\n"
               "    return it",
      py="mandelgrid(100)", py_check="int(mandelgrid(100).sum())")
bench(id="mcpi", group=SC, name="Monte Carlo pi, 10^7 samples (vectorized)",
      setup="n=10^7;",
      expr="4. Total[UnitStep[1. - RandomReal[{0,1},n]^2 "
           "- RandomReal[{0,1},n]^2]]/n;", reps=3,
      py_setup="n=10**7",
      py="4.0*np.count_nonzero(1.0 - np.random.rand(n)**2 "
         "- np.random.rand(n)**2 >= 0.0)/n")
# The point set is a cubic LATTICE, deliberately. A pseudo-random one built from
# Mod[a k, m] cycles, so with 3n coordinates it repeats POINTS; two coincident
# bodies make r = 0, the 1/r^12 term overflows to Inf, and the compiled call then
# correctly bails to the interpreter on the non-finite result -- which read as an
# 89x performance cliff between n=1000 and n=1200 and answered Indeterminate.
# A lattice with spacing 1.2 has no coincident pairs and a finite energy.
bench(id="lj", group=SC, name="Lennard-Jones energy, 1452 bodies (all pairs)",
      setup="nb=1452; pts=Flatten[Table[1.2 {i,j,k2}, "
            "{i,0,10},{j,0,10},{k2,0,11}]]; "
            "lj=Compile[{{p,_Real,1},{m,_Integer}}, "
            "Module[{e=0.,i=0,j=0,dx=0.,dy=0.,dz=0.,r2=0.,ir6=0.}, "
            "For[i=0,i<m-1,i++, For[j=i+1,j<m,j++, "
            "dx=p[[3i+1]]-p[[3j+1]]; dy=p[[3i+2]]-p[[3j+2]]; "
            "dz=p[[3i+3]]-p[[3j+3]]; r2=dx dx+dy dy+dz dz; "
            "ir6=1./(r2 r2 r2); e=e+4. ir6 (ir6-1.)]]; e]];",
      expr="lj[pts, nb];", check="lj[pts, nb]", reps=3,
      py_setup="nb=1452\n"
               "pts=np.array([1.2*c for i in range(11) for j in range(11)\n"
               "              for k in range(12) for c in (i,j,k)])\n"
               "iu=np.triu_indices(nb,1)\n"
               "def lj(p, m):\n"
               "    q = p.reshape(m,3)\n"
               "    d = q[:,None,:] - q[None,:,:]\n"
               "    r2 = (d*d).sum(-1)[iu]\n"
               "    ir6 = 1.0/(r2*r2*r2)\n"
               "    return float((4.0*ir6*(ir6-1.0)).sum())",
      py="lj(pts, nb)", py_check="lj(pts, nb)")

# ---- integer / combinatorial --------------------------------------------
bench(id="sieve", group=IN, name="Sieve of Eratosthenes to 10^7",
      setup="sv=Compile[{{m,_Integer}}, Module[{s=ConstantArray[1,m],i=2,j=0,c=0}, "
            "While[i i<=m, If[s[[i]]==1, j=i i; While[j<=m, s[[j]]=0; j=j+i]]; "
            "i=i+1]; i=2; While[i<=m, c=c+s[[i]]; i=i+1]; c]];",
      expr="sv[10^7];", check="sv[10^7]", reps=3,
      py_setup="def sv(m):\n"
               "    s = np.ones(m+1, dtype=bool); s[:2] = False\n"
               "    i = 2\n"
               "    while i*i <= m:\n"
               "        if s[i]: s[i*i::i] = False\n"
               "        i += 1\n"
               "    return int(s.sum())",
      py="sv(10**7)", py_check="sv(10**7)")
bench(id="collatz", group=IN, name="Collatz longest chain below 10^6",
      setup="cz=Compile[{{m,_Integer}}, Module[{bl=0,k=1,q=0,len=0}, "
            "While[k<=m, q=k; len=0; While[q>1, "
            "If[Mod[q,2]==0, q=Quotient[q,2], q=3q+1]; len=len+1]; "
            "If[len>bl, bl=len]; k=k+1]; bl]];",
      expr="cz[10^6];", check="cz[10^6]", reps=3,
      py_setup="def cz(m):\n"
               "    bl = 0\n"
               "    for k in range(1, m+1):\n"
               "        q = k; ln = 0\n"
               "        while q > 1:\n"
               "            q = q//2 if q % 2 == 0 else 3*q+1\n"
               "            ln += 1\n"
               "        if ln > bl: bl = ln\n"
               "    return bl",
      py="cz(10**6)", py_check="cz(10**6)",
      note="Python column is CPython, not NumPy: no numba on this host")
bench(id="fib", group=IN, name="Naive recursive Fibonacci, fib(25)",
      setup="Clear[fib]; fib[0]=0; fib[1]=1; fib[k_]:=fib[k-1]+fib[k-2];",
      expr="fib[25];", check="fib[25]", reps=3,
      py_setup="def fib(k):\n    return k if k < 2 else fib(k-1)+fib(k-2)",
      py="fib(25)", py_check="fib(25)",
      note="rule dispatch, not compiled in any of the three")
bench(id="primepi", group=IN, name="PrimePi[10^9]", setup="",
      expr="PrimePi[10^9 + kk];", check="PrimePi[10^9]", reps=3)

# ---- arbitrary precision ------------------------------------------------
# NumPy has no arbitrary precision; the Python column uses mpmath (which is what
# a Python user reaching for 100,000 digits would use) and CPython's own bignums,
# and is labelled accordingly.
bench(id="pi", group=BN, name="pi to 100,000 digits", setup="",
      expr="N[Pi,100000 + kk];", check="Round[10^20 (N[Pi,100000] - 3)]", reps=3,
      cold=True,
      py_setup="import mpmath as mp\n"
               "def mppi(d):\n"
               "    mp.mp.dps = d + 10\n"
               "    return +mp.pi",
      py="mppi(100000 + kk)",
      py_check="int(mp.nint((mppi(100000) - 3) * mp.mpf(10)**20))",
      note="Python column is mpmath, not NumPy")
bench(id="fact", group=BN, name="50000! (exact)", setup="",
      expr="(50000 + kk)!;", check="IntegerLength[50000!]", reps=3,
      # set_int_max_str_digits here as well as in bigmul: CPython caps int->str
      # at 4300 digits by default, this row runs FIRST, and the cap made its
      # check raise while bigmul's identical check passed.
      py_setup="import math as _m\nimport sys as _s\n_s.set_int_max_str_digits(0)",
      py="_m.factorial(50000 + kk)",
      py_check="len(str(_m.factorial(50000)))",
      note="Python column is CPython bignums, not NumPy")
# IntegerLength[], not `expr;`. Wolfram does not materialise a discarded bignum
# product: `(p + kk) q;` reads 138 us and the forced form 2.96 ms, which would
# have been reported as a 24x gap that does not exist. Array ops were checked in
# both forms and agree within noise; bignum arithmetic is the exception.
#
# The Python side has the mirror-image trap: len(str(n)) is the natural spelling
# of IntegerLength but CPython's int->decimal conversion is quadratic, so timing
# it would measure base conversion rather than the multiply. The timed Python
# expression forces the product with .bit_length() (O(1) on a materialised int)
# and only the untimed check pays for the decimal form.
bench(id="bigmul", group=BN, name="Product of two 10^6-bit integers",
      setup="p=2^1000003 - 1; q=3^631305;",
      expr="IntegerLength[(p + kk) q];", check="IntegerLength[p q]", reps=5,
      py_setup="import sys as _s\n_s.set_int_max_str_digits(0)\n"
               "p=2**1000003 - 1\nq=3**631305",
      py="((p + kk)*q).bit_length()", py_check="len(str(p*q))",
      note="Python column is CPython bignums, not NumPy")
bench(id="bigpow", group=BN, name="3^1000000", setup="",
      expr="3^(1000000 + kk);", check="IntegerLength[3^1000000]", reps=3,
      py_setup="import sys as _s\n_s.set_int_max_str_digits(0)",
      py="3**(1000000 + kk)", py_check="len(str(3**1000000))",
      note="Python column is CPython bignums, not NumPy")

# --------------------------------------------------------------------------
# Second sweep (2026-07-31). Five kernels chosen to probe subsystems the first
# 38 never touched: a Krylov iterative solver, direct convolution, ODE time
# integration, an irregular hash-keyed reduction, and interpolation. The point
# of a new benchmark is to run somewhere the old ones did not.
# --------------------------------------------------------------------------
IT = "Iterative solvers and convolution"
OD = "ODE integration and interpolation"
IR = "Irregular / data-dependent (n = 10^7)"

# Matrix-free conjugate gradient on the periodic 2D Poisson operator: the HPCG
# shape, and the composition test the Jacobi row is not -- a Krylov iteration
# alternates stencil applications with GLOBAL reductions (Total[p ap, 2]), so a
# single unpacked temporary anywhere in the loop shows up.
#
# The source is a POINT SOURCE with the mean removed, and both halves of that
# matter. The periodic Laplacian is singular with the constants in its
# nullspace, so a source with a nonzero mean makes the solve ill-posed; and a
# smooth source such as Sin[2 Pi i/n] Cos[4 Pi j/n] is an EIGENVECTOR of the
# discrete operator, which CG annihilates in a single iteration -- the residual
# reached 3.6e-24 in 10 iterations and the remaining 90 measured denormal
# arithmetic. A delta has a broad spectrum, so 100 iterations is 100 iterations
# of genuine work and the residual is a meaningful cross-system check.
bench(id="cg", group=IT, name="Conjugate gradient, 2D Poisson 256^2, 100 iterations",
      # The delta is placed and the mean removed IN TERMS OF n, so the problem
      # stays well-posed under --scale. Written with the literals 128 and 65536
      # it did not: at a shrunken n the If never fired, b became a constant with
      # a nonzero mean, and the periodic operator is singular on the constants --
      # so Total[p ap, 2] was 0, alpha was ComplexInfinity, and the benchmark
      # measured Indeterminate propagating through an NDArray.
      setup="n=256; c0=Floor[n/2]; b=Table[If[i==c0 && j==c0, 1.-1./n^2, -1./n^2],{i,n},{j,n}]; "
            "lap[u_]:=4. u - RotateLeft[u,{1,0}] - RotateRight[u,{1,0}] "
            "- RotateLeft[u,{0,1}] - RotateRight[u,{0,1}]; "
            "cg[m_]:=Module[{x,r,p,ap,rs,rsn,al,k}, x=ConstantArray[0.,{n,n}]; "
            "r=b; p=r; rs=Total[r r,2]; k=0; While[k<m, ap=lap[p]; "
            "al=rs/Total[p ap,2]; x=x+al p; r=r-al ap; rsn=Total[r r,2]; "
            "p=r+(rsn/rs) p; rs=rsn; k=k+1]; rs];",
      expr="cg[100];", check="cg[100]", reps=3,
      py_setup="n=256\nc0=n//2\n"
               "b=np.full((n,n), -1.0/n**2); b[c0-1,c0-1] = 1.0-1.0/n**2\n"
               "def lap(u):\n"
               "    return (4.0*u - np.roll(u,-1,0) - np.roll(u,1,0)\n"
               "            - np.roll(u,-1,1) - np.roll(u,1,1))\n"
               "def cg(m):\n"
               "    x = np.zeros((n,n)); r = b.copy(); p = r.copy()\n"
               "    rs = (r*r).sum()\n"
               "    for _ in range(m):\n"
               "        ap = lap(p)\n"
               "        al = rs/(p*ap).sum()\n"
               "        x = x + al*p; r = r - al*ap\n"
               "        rsn = (r*r).sum(); p = r + (rsn/rs)*p; rs = rsn\n"
               "    return float(rs)",
      py="cg(100)", py_check="cg(100)")

# Direct 2D convolution -- image filtering, and the same arithmetic as a stencil
# written the other way round. The timed data is random; the CHECK runs a small
# deterministic instance instead, because the two systems cannot be made to draw
# the same random image and a value check is worth more than a matching one.
bench(id="conv2d", group=IT, name="ListConvolve, 1024^2 image, 5x5 kernel",
      setup="img=RandomReal[{0,1},{1024,1024}]; ker=RandomReal[{0,1},{5,5}];",
      expr="ListConvolve[ker,img];",
      check="Total[ListConvolve[{{1.,2.},{3.,4.}}, "
            "Table[N[Mod[i j,7]],{i,64},{j,64}]],2]", reps=3,
      py_setup="img=np.random.rand(1024,1024)\nker=np.random.rand(5,5)\n"
               "cimg=np.array([[float((i*j) % 7) for j in range(1,65)]\n"
               "               for i in range(1,65)])",
      py="spsig.convolve(img, ker, mode='valid')",
      py_check="float(spsig.convolve(cimg, np.array([[1.,2.],[3.,4.]]), "
               "mode='valid').sum())")

# Lorenz to t = 200: adaptive-step ODE integration, the workload NDSolve exists
# for. The CHECK is a DIFFERENT problem on purpose -- Lorenz is chaotic, so its
# state at t = 200 is a property of the step sequence and not of the equations,
# and requiring two independently-written solvers to agree on it would be
# meaningless. y' = -y on [0,1] has the answer 1/E, all three solvers target ~8
# digits, and it verifies the same code path end to end.
#
# The dependent variables are lx/ly/lz, NOT x/y/z, and that is load-bearing: the
# array section above binds `x` and `y` to 10^7-element lists, every benchmark
# here shares one kernel session, and an ODE whose dependent variable already has
# a value is a different problem. Wolfram solved it anyway and read 4x slower;
# Mathilda SEGFAULTED. (The crash is not caused by this file -- it reproduces on
# a clean checkout from `x = RandomReal[{0,1}, 100000]` followed by an NDSolve in
# `x`, and is recorded in the HPC plan. The benchmark should not have been
# stepping on a bound symbol regardless.)
bench(id="ode", group=OD, name="NDSolve, Lorenz system to t = 200",
      setup="Clear[lx,ly,lz];",
      expr="NDSolve[{lx'[t]==10.(ly[t]-lx[t]),ly'[t]==lx[t](28.-lz[t])-ly[t],"
           "lz'[t]==lx[t] ly[t]-(8./3.) lz[t],lx[0]==1.,ly[0]==1.,lz[0]==1.},"
           "{lx,ly,lz},{t,0,200}];",
      check="yy[1] /. First[NDSolve[{yy'[t]==-yy[t],yy[0]==1.},yy,{t,0,1}]]",
      reps=3,
      py_setup="from scipy.integrate import solve_ivp\n"
               "def lorenz(t, s):\n"
               "    return [10.0*(s[1]-s[0]), s[0]*(28.0-s[2])-s[1],\n"
               "            s[0]*s[1]-(8.0/3.0)*s[2]]",
      py="solve_ivp(lorenz, (0.0,200.0), [1.0,1.0,1.0], rtol=1e-8, atol=1e-8)",
      py_check="float(solve_ivp(lambda t,s: -s, (0.0,1.0), [1.0], "
               "rtol=1e-9, atol=1e-9).y[0][-1])")

# Build a cubic interpolation and evaluate it at as many points as it has nodes
# -- resampling/regridding, and the standard way a tabulated function is used.
# Deliberately modest sizes: the per-evaluation cost turned out to scale with
# the TABLE size in Mathilda, so 10^6 nodes could not be timed at all.
#
# No py_check: SciPy's CubicSpline is a global C2 spline and Wolfram's
# Interpolation is a local piecewise fit, so the two do not agree at a point and
# should not be asked to. The timings are still comparable -- both build a cubic
# interpolant over the same nodes and evaluate it at the same count -- and the
# Mathilda/Wolfram value check is unaffected.
bench(id="interp", group=OD, name="Interpolation, 10^4 nodes, 10^4 evaluations",
      setup="ni=10^4; tbl=Table[{xx,Sin[xx]},{xx,0.,100.,100./(ni-1)}]; "
            "ifn=Interpolation[tbl]; pts=RandomReal[{0,100},ni];",
      expr="ifn /@ pts;", check="ifn[3.3]", reps=3,
      py_setup="from scipy.interpolate import CubicSpline\n"
               "ni=10**4\nxs=np.linspace(0.0,100.0,ni)\n"
               "ifn=CubicSpline(xs, np.sin(xs))\n"
               "pts=np.random.rand(ni)*100.0",
      py="ifn(pts)")

# Tally: an irregular, data-dependent reduction -- a scatter-add into a keyed
# table, where the streaming kernels of section 4 all have a static access
# pattern. With 10^7 draws over 10^4 bins every bin is hit with probability
# indistinguishable from 1, so Length[] is a deterministic check despite the
# input being random.
bench(id="tally", group=IR, name="Tally, 10^7 integers into 10^4 bins",
      setup="nt=10^7; ii=RandomInteger[{1,10^4},nt];",
      expr="Tally[ii];", check="Length[Tally[ii]]", reps=3,
      py_setup="nt=10**7\nii=np.random.randint(1,10**4+1,nt)",
      py="np.unique(ii, return_counts=True)",
      py_check="int(np.unique(ii).size)")

# --------------------------------------------------------------------------
# Third sweep (2026-07-31). Eight kernels chosen for DEPTH rather than breadth:
# each is a pipeline out of a real application, not a single builtin. The second
# sweep's CG row is the reason -- a composition can be slow in ways that none of
# its parts are, and only a composition finds that.
#
# Every one runs in all three systems, and every one has a deterministic scalar
# check that all three agree on. That is a stronger constraint than it sounds:
# the check pins the ALGORITHM, so the three columns cannot silently be timing
# three different computations.
# --------------------------------------------------------------------------
FI = "Computational finance"
ML = "Machine learning"
NB = "N-body and CFD"
SG = "Signal and image processing"

# Black-Scholes by Monte Carlo. The whole of the option-pricing workload: draw,
# transform to normals (Box-Muller: Log, Sqrt, Cos over 10^7), map through the
# lognormal, take a payoff, reduce.
#
# The check is the CLOSED FORM, not the Monte Carlo estimate -- a random estimate
# cannot be compared across systems that do not share a generator, and the
# analytic price is both deterministic and a real test of Erf. At s0 = strike =
# 100, r = 0.05, sigma = 0.2, T = 1 the d's are exactly 0.35 and 0.15, so the
# check needs no intermediate algebra and can be read at a glance.
bench(id="bsmc", group=FI, name="Black-Scholes Monte Carlo, 10^7 paths",
      setup="npth=10^7; s0=100.; strike=100.; rr=0.05; sig=0.2; tt=1.;"
            "bs[m_]:=Module[{u1,u2,z,st}, u1=RandomReal[{0,1},m]; "
            "u2=RandomReal[{0,1},m]; "
            "z=Sqrt[-2. Log[u1]] Cos[6.283185307179586 u2]; "
            "st=s0 Exp[(rr-0.5 sig^2) tt + sig Sqrt[tt] z]; "
            "Exp[-rr tt] Total[UnitStep[st-strike](st-strike)]/m];",
      expr="bs[npth];",
      check="100. 0.5 (1.+Erf[0.35/Sqrt[2.]]) "
            "- 100. Exp[-0.05] 0.5 (1.+Erf[0.15/Sqrt[2.]])",
      reps=3,
      py_setup="npth=10**7\ns0=100.0\nstrike=100.0\nrr=0.05\nsig=0.2\ntt=1.0\n"
               "def bs(m):\n"
               "    u1 = np.random.rand(m); u2 = np.random.rand(m)\n"
               "    z = np.sqrt(-2.0*np.log(u1))*np.cos(6.283185307179586*u2)\n"
               "    st = s0*np.exp((rr-0.5*sig**2)*tt + sig*np.sqrt(tt)*z)\n"
               "    pay = np.where(st-strike >= 0.0, st-strike, 0.0)\n"
               "    return float(math.exp(-rr*tt)*pay.sum()/m)",
      py="bs(npth)",
      py_check="100.0*0.5*(1.0+math.erf(0.35/math.sqrt(2.0))) "
               "- 100.0*math.exp(-0.05)*0.5*(1.0+math.erf(0.15/math.sqrt(2.0)))")

# A quant time-series pipeline: cumulative log return -> price path, an
# exponential moving average, a rolling volatility, and the maximum drawdown.
#
# This one is here because two of the four stages are SCANS, and a scan is the
# thing an array language cannot vectorize. Accumulate covers the prefix sum, but
# the EMA is a general linear recurrence and the drawdown needs a running
# maximum, and neither has a builtin in either CAS -- both are written as
# FoldList, i.e. one interpreted call per element. NumPy is in the same position
# and its answer is scipy.signal.lfilter and np.maximum.accumulate, which is what
# the row measures: the cost of a scan done in C versus a scan done in an
# evaluator. The rolling volatility is deliberately written the vectorizable way
# (differences of a prefix sum) so the row is not purely about the scans.
bench(id="tseries", group=FI, name="Return series: EMA + rolling vol + drawdown, 10^6",
      setup="nts=10^6; wts=250; ret=RandomReal[{-0.01,0.01},nts]; "
            "tsr[r_,w_]:=Module[{px,ema,run,cs,vol}, px=Exp[Accumulate[r]]; "
            "ema=FoldList[0.98 #1+0.02 #2&,First[px],Rest[px]]; "
            "run=FoldList[Max,First[px],Rest[px]]; cs=Accumulate[r r]; "
            "vol=Sqrt[(Drop[cs,w]-Drop[cs,-w])/w]; "
            "Last[ema]+Max[1.-px/run]+Max[vol]];",
      expr="tsr[ret,wts];",
      check="tsr[Table[N[Mod[3 i, 7]]/700. - 0.004, {i,2000}], 50]", reps=3,
      py_setup="from scipy.signal import lfilter\n"
               "nts=10**6\nwts=250\nret=np.random.rand(nts)*0.02-0.01\n"
               "def tsr(r, w):\n"
               "    px = np.exp(np.cumsum(r))\n"
               "    ema = lfilter([0.02],[1.0,-0.98], px[1:], zi=[0.98*px[0]])[0]\n"
               "    run = np.maximum.accumulate(px)\n"
               "    cs = np.cumsum(r*r)\n"
               "    vol = np.sqrt((cs[w:]-cs[:-w])/w)\n"
               "    return (float(ema[-1]) + float((1.0-px/run).max())\n"
               "            + float(vol.max()))",
      py="tsr(ret,wts)",
      py_check="tsr(np.array([float((3*i) % 7)/700.0-0.004 "
               "for i in range(1,2001)]), 50)")

# Logistic regression by batch gradient descent: the inner loop of essentially
# every linear model. Two matrix-vector products per iteration in OPPOSITE
# orientations (X.w and Transpose[X].e), a transcendental over the full sample,
# and loop-carried state -- so the temporaries cannot be hoisted and a single
# unpacked intermediate is paid 100 times.
bench(id="logreg", group=ML, name="Logistic regression, 200000x32, 100 GD steps",
      setup="nlr=200000; dlr=32; Xlr=RandomReal[{-1,1},{nlr,dlr}]; "
            "wtr=RandomReal[{-1,1},dlr]; ylr=UnitStep[Xlr . wtr]; "
            "lrfit[xm_,yv_,m_]:=Module[{w,pr,g,k,nn,dd}, nn=Length[xm]; "
            "dd=Length[xm[[1]]]; w=ConstantArray[0.,dd]; k=0; "
            "While[k<m, pr=1./(1.+Exp[-(xm . w)]); "
            "g=(Transpose[xm] . (pr-yv))/nn; w=w-4. g; k=k+1]; Total[w w]];",
      expr="lrfit[Xlr,ylr,100];",
      check="lrfit[Table[N[Mod[i k, 5]]/5. - 0.4, {i,64},{k,4}], "
            "Table[N[Mod[i,2]], {i,64}], 50]", reps=3,
      py_setup="nlr=200000\ndlr=32\nXlr=np.random.rand(nlr,dlr)*2.0-1.0\n"
               "wtr=np.random.rand(dlr)*2.0-1.0\n"
               "ylr=(Xlr @ wtr >= 0.0).astype(float)\n"
               "def lrfit(xm, yv, m):\n"
               "    nn, dd = xm.shape\n"
               "    w = np.zeros(dd)\n"
               "    for _ in range(m):\n"
               "        pr = 1.0/(1.0+np.exp(-(xm @ w)))\n"
               "        w = w - 4.0*((xm.T @ (pr-yv))/nn)\n"
               "    return float(w @ w)\n"
               "Xchk=np.array([[float((i*k) % 5)/5.0-0.4 for k in range(1,5)]\n"
               "               for i in range(1,65)])\n"
               "ychk=np.array([float(i % 2) for i in range(1,65)])",
      py="lrfit(Xlr,ylr,100)", py_check="lrfit(Xchk,ychk,50)")

# Lloyd's algorithm, written the way an array language has to write it: the
# assignment step is a k x n distance matrix and an elementwise minimum across
# it, and the update step is a DOT PRODUCT against a 0/1 mask -- which is how a
# scatter-add is spelled without a scatter.
#
# The masks are made mutually exclusive by subtracting what has already been
# assigned, rather than by trusting that no two centroids tie. With ties allowed
# a point would be counted in two clusters and the result would depend on
# floating-point luck; the check is a deterministic dataset, so that would not
# have been a theoretical concern.
bench(id="kmeans", group=ML, name="k-means, 100000 points, 8 dims, k = 16, 20 its",
      setup="nkm=100000; dkm=8; kkm=16; "
            "Pkm=Transpose[RandomReal[{0,1},{nkm,dkm}]]; "
            "km[ptsT_,kc_,m_]:=Module[{cs,d,mnv,asg,mj,cnt,i,nn}, "
            "nn=Length[ptsT[[1]]]; cs=Table[ptsT[[All,j]],{j,kc}]; i=0; "
            "While[i<m, d=Table[Total[(ptsT-cs[[j]])^2],{j,kc}]; "
            "mnv=MapThread[Min,d]; asg=ConstantArray[0.,nn]; "
            "cs=Table[mj=UnitStep[mnv-d[[j]]](1.-asg); asg=asg+mj; "
            "cnt=Total[mj]; If[cnt>0.,(ptsT . mj)/cnt,cs[[j]]],{j,kc}]; "
            "i=i+1]; "
            "Total[MapThread[Min,Table[Total[(ptsT-cs[[j]])^2],{j,kc}]]]];",
      expr="km[Pkm,kkm,20];",
      check="km[Transpose[Table[N[Mod[i^2+3 i k, 101]]/101.+0.001 N[k],"
            "{i,200},{k,4}]], 4, 10]", reps=3,
      py_setup="nkm=100000\ndkm=8\nkkm=16\n"
               "Pkm=np.random.rand(nkm,dkm).T.copy()\n"
               "def km(ptsT, kc, m):\n"
               "    nn = ptsT.shape[1]\n"
               "    cs = [ptsT[:,j].copy() for j in range(kc)]\n"
               "    def dist(c):\n"
               "        return np.array([((ptsT-q[:,None])**2).sum(axis=0)"
               " for q in c])\n"
               "    for _ in range(m):\n"
               "        d = dist(cs); mnv = d.min(axis=0)\n"
               "        asg = np.zeros(nn); new = []\n"
               "        for j in range(kc):\n"
               "            mj = np.where(mnv-d[j] >= 0.0, 1.0, 0.0)*(1.0-asg)\n"
               "            asg = asg + mj; cnt = mj.sum()\n"
               "            new.append((ptsT @ mj)/cnt if cnt > 0.0 else cs[j])\n"
               "        cs = new\n"
               "    return float(dist(cs).min(axis=0).sum())\n"
               "Pchk=np.array([[float((i*i+3*i*k) % 101)/101.0+0.001*k\n"
               "                for k in range(1,5)] for i in range(1,201)]).T.copy()",
      py="km(Pkm,kkm,20)", py_check="km(Pchk,4,10)")

# All-pairs gravity with a velocity-Verlet step, written in the ARRAY style: the
# pairwise displacements are three n x n matrices from Outer[Subtract, .., ..].
# This is here specifically for Outer, which was noticed at ~440 ns/element while
# looking at something else and never benchmarked. 1024 bodies makes each of
# those matrices a megabyte-scale intermediate, six of them per step.
bench(id="nbody", group=NB, name="N-body all-pairs gravity, 1024 bodies, 10 steps",
      setup="nbd=1024; "
            "st0=Table[RandomReal[{0,1},nbd],{6}]; "
            "grav[s_]:=Module[{px,py,pz,vx,vy,vz,dx,dy,dz,r2,iv,ax,ay,az,"
            "dt=0.001,eps=0.01}, px=s[[1]];py=s[[2]];pz=s[[3]]; "
            "vx=s[[4]];vy=s[[5]];vz=s[[6]]; "
            "dx=Outer[Subtract,px,px]; dy=Outer[Subtract,py,py]; "
            "dz=Outer[Subtract,pz,pz]; r2=dx dx+dy dy+dz dz+eps; "
            "iv=r2^(-1.5); ax=-Total[dx iv,{2}]; ay=-Total[dy iv,{2}]; "
            "az=-Total[dz iv,{2}]; "
            "{px+dt(vx+0.5 dt ax),py+dt(vy+0.5 dt ay),pz+dt(vz+0.5 dt az),"
            "vx+dt ax,vy+dt ay,vz+dt az}];",
      expr="Nest[grav,st0,10];",
      check="Total[Nest[grav, Table[N[Mod[7 i + 3 k, 11]]/11., {i,6},{k,8}], 3], 2]",
      reps=3,
      py_setup="nbd=1024\nst0=[np.random.rand(nbd) for _ in range(6)]\n"
               "def grav(s, dt=0.001, eps=0.01):\n"
               "    px,py,pz,vx,vy,vz = s\n"
               "    dx = px[:,None]-px[None,:]\n"
               "    dy = py[:,None]-py[None,:]\n"
               "    dz = pz[:,None]-pz[None,:]\n"
               "    r2 = dx*dx+dy*dy+dz*dz+eps\n"
               "    iv = r2**-1.5\n"
               "    ax = -(dx*iv).sum(axis=1); ay = -(dy*iv).sum(axis=1)\n"
               "    az = -(dz*iv).sum(axis=1)\n"
               "    return [px+dt*(vx+0.5*dt*ax), py+dt*(vy+0.5*dt*ay),\n"
               "            pz+dt*(vz+0.5*dt*az), vx+dt*ax, vy+dt*ay, vz+dt*az]\n"
               "schk=[np.array([float((7*i+3*k) % 11)/11.0 for k in range(1,9)])\n"
               "      for i in range(1,7)]",
      py="nest(grav,st0,10)",
      py_check="float(sum(a.sum() for a in nest(grav,schk,3)))")

# A 3D 7-point heat stencil. The stencil rows above are all rank 2; this one is
# rank 3, where the shift is along an axis whose stride is neither 1 nor the row
# length, and the working set (128^3 doubles = 16 MB) exceeds any cache on this
# machine. Six shifted copies per step, 50 steps.
bench(id="heat3d", group=NB, name="3D heat equation, 128^3, 50 steps",
      setup="n3=128; u3=RandomReal[{0,1},{n3,n3,n3}]; "
            "hstep[u_]:=u+0.1(RotateLeft[u,{1,0,0}]+RotateRight[u,{1,0,0}]"
            "+RotateLeft[u,{0,1,0}]+RotateRight[u,{0,1,0}]"
            "+RotateLeft[u,{0,0,1}]+RotateRight[u,{0,0,1}]-6. u);",
      expr="Nest[hstep,u3,50];",
      check="Total[Nest[hstep, Table[N[Mod[i+2j+3k,5]],{i,8},{j,8},{k,8}], 5], 3]",
      reps=3,
      py_setup="n3=128\nu3=np.random.rand(n3,n3,n3)\n"
               "def hstep(u):\n"
               "    return u + 0.1*(np.roll(u,-1,0)+np.roll(u,1,0)\n"
               "                    +np.roll(u,-1,1)+np.roll(u,1,1)\n"
               "                    +np.roll(u,-1,2)+np.roll(u,1,2)-6.0*u)\n"
               "uchk=np.array([[[float((i+2*j+3*k) % 5) for k in range(1,9)]\n"
               "                for j in range(1,9)] for i in range(1,9)])",
      py="nest(hstep,u3,50)",
      py_check="float(nest(hstep,uchk,5).sum())")

# Welch's method: split a signal into blocks, window each, transform, accumulate
# the squared magnitude. The standard way a power spectrum is actually computed,
# and a different shape of FFT work from the two single-transform rows above --
# 1024 medium transforms with windowing and accumulation between them, where
# those are one enormous transform in isolation.
#
# The explicit /Sqrt[bl] is what makes the three columns comparable: Fourier[]
# carries a 1/Sqrt[n] and numpy.fft.fft carries none, so without it the check
# would differ by a factor of the block size and the row would be measuring the
# same work while reporting different answers.
bench(id="psd", group=SG, name="Welch PSD, 2^22 samples, 1024 blocks of 4096",
      setup="nsg=2^22; blk=4096; nblk=1024; "
            "sgn=RandomReal[{-1,1},nsg]; "
            "wnd=Table[0.5(1.-Cos[2. Pi (j-1)/(blk-1)]),{j,blk}]; "
            "psd[s_,w_,nb_,bl_]:=Module[{acc,i,seg}, acc=ConstantArray[0.,bl]; "
            "i=0; While[i<nb, seg=Take[s,{i bl+1,i bl+bl}] w; "
            "acc=acc+Abs[Fourier[seg]]^2; i=i+1]; acc];",
      expr="psd[sgn,wnd,nblk,blk];",
      check="Total[psd[Table[N[Sin[0.1 j]],{j,4096}], "
            "Table[0.5(1.-Cos[2. Pi (j-1)/255.]),{j,256}], 16, 256]]",
      reps=3,
      py_setup="nsg=2**22\nblk=4096\nnblk=1024\n"
               "sgn=np.random.rand(nsg)*2.0-1.0\n"
               "wnd=0.5*(1.0-np.cos(2.0*math.pi*np.arange(blk)/(blk-1)))\n"
               "def psd(s, w, nb, bl):\n"
               "    acc = np.zeros(bl)\n"
               "    rt = math.sqrt(bl)\n"
               "    for i in range(nb):\n"
               "        acc = acc + np.abs(np.fft.fft(s[i*bl:i*bl+bl]*w)/rt)**2\n"
               "    return acc\n"
               "schk=np.sin(0.1*np.arange(1,4097))\n"
               "wchk=0.5*(1.0-np.cos(2.0*math.pi*np.arange(0,256)/255.0))",
      py="psd(sgn,wnd,nblk,blk)",
      py_check="float(psd(schk,wchk,16,256).sum())")

# An edge-detection pipeline: separable 5-tap Gaussian blur (two rank-2
# correlations with a 1x5 and a 5x1 kernel), then the two Sobel derivatives, then
# the gradient magnitude. Four correlations and a Sqrt over a megapixel.
#
# The kernels are 1x5 and 5x1 rather than length-5 vectors because
# ListCorrelate[rank1, rank2] is unevaluated in BOTH Mathilda and Wolfram
# (ListCorrelate::kldims) -- checked, not assumed, and not a Mathilda gap.
bench(id="imgpipe", group=SG, name="Gaussian blur + Sobel edges, 1024^2",
      setup="imr=RandomReal[{0,1},{1024,1024}]; "
            "gk={{1.,4.,6.,4.,1.}}/16.; gkT=Transpose[gk]; "
            "sxk={{-1.,0.,1.},{-2.,0.,2.},{-1.,0.,1.}}; "
            "syk={{-1.,-2.,-1.},{0.,0.,0.},{1.,2.,1.}}; "
            "edge[im_]:=Module[{bl,gx,gy}, "
            "bl=ListCorrelate[gkT,ListCorrelate[gk,im]]; "
            "gx=ListCorrelate[sxk,bl]; gy=ListCorrelate[syk,bl]; "
            "Sqrt[gx gx+gy gy]];",
      expr="edge[imr];",
      check="Total[edge[Table[N[Mod[i j, 17]],{i,64},{j,64}]], 2]", reps=3,
      py_setup="imr=np.random.rand(1024,1024)\n"
               "gk=np.array([[1.,4.,6.,4.,1.]])/16.0\ngkT=gk.T.copy()\n"
               "sxk=np.array([[-1.,0.,1.],[-2.,0.,2.],[-1.,0.,1.]])\n"
               "syk=np.array([[-1.,-2.,-1.],[0.,0.,0.],[1.,2.,1.]])\n"
               "def edge(im):\n"
               "    b = spsig.correlate(spsig.correlate(im, gk, mode='valid'),\n"
               "                      gkT, mode='valid')\n"
               "    gx = spsig.correlate(b, sxk, mode='valid')\n"
               "    gy = spsig.correlate(b, syk, mode='valid')\n"
               "    return np.sqrt(gx*gx+gy*gy)\n"
               "ichk=np.array([[float((i*j) % 17) for j in range(1,65)]\n"
               "               for i in range(1,65)])",
      py="edge(imr)", py_check="float(edge(ichk).sum())")

SCALED = {  # id -> (pattern, replacement-fn) applied to setup+expr under --scale
    "matmul": [("n=1000", "n")], "solve": [("n=1000", "n")],
    "inverse": [("n=500", "n")], "det": [("n=500", "n")],
    "qr": [("n=500", "n")], "eigen": [("n=300", "n")], "svd": [("n=300", "n")],
}


def apply_scale(b, s):
    """Shrink the dominant size knob so the whole suite can smoke-test fast."""
    if s == 1.0:
        return b
    b = dict(b)
    for key in ("setup", "expr", "check", "wl", "wl_setup",
                "py", "py_setup", "py_check"):
        t = b.get(key)
        if not t:
            continue
        def shrink(m):
            v = int(m.group(2))
            return m.group(1) + str(max(2, int(v * (s ** 0.5) if v > 64 else v)))
        t = re.sub(r"(n=|nb=|nbd=|n3=|m=)(\d+)", shrink, t)
        # Both spellings of every power: the CAS say 10^7, Python says 10**7.
        for a, c in (("10^7", "10^5"), ("10**7", "10**5"),
                     ("10^6", "10^4"), ("10**6", "10**4"),
                     ("2^20", "2^14"), ("2**20", "2**14"),
                     ("2^22", "2^16"), ("2**22", "2**16"),
                     ("10^9", "10^6"), ("10**9", "10**6"),
                     ("3^1000000", "3^20000"), ("3**1000000", "3**20000")):
            t = t.replace(a, c)
        t = t.replace("100000", "2000").replace("50000", "2000")
        t = t.replace(",100]", ",20]")
        t = t.replace("u0,100]", "u0,5]").replace("g0,100]", "g0,5]")
        t = t.replace("u0,100)", "u0,5)").replace("g0,100)", "g0,5)")
        t = t.replace("cg[100]", "cg[5]").replace("cg(100)", "cg(5)")
        t = t.replace("{t,0,200}", "{t,0,5}").replace("(0.0,200.0)", "(0.0,5.0)")
        t = t.replace("{1024,1024}", "{128,128}").replace("(1024,1024)", "(128,128)")
        t = t.replace("ni=10^4", "ni=500").replace("ni=10**4", "ni=500")
        t = t.replace("nblk=1024", "nblk=32")
        t = t.replace("u3,50]", "u3,5]").replace("u3,50)", "u3,5)")
        t = t.replace("st0,10]", "st0,3]").replace("st0,10)", "st0,3)")
        t = t.replace("Xlr,ylr,100]", "Xlr,ylr,10]")
        t = t.replace("Xlr,ylr,100)", "Xlr,ylr,10)")
        t = t.replace("Pkm,kkm,20]", "Pkm,kkm,3]").replace("Pkm,kkm,20)", "Pkm,kkm,3)")
        b[key] = t
    return b


# --------------------------------------------------------------------------
# Mathilda: NDJSON pipe protocol (non-TTY stdin).
# --------------------------------------------------------------------------
def run_mathilda(benches):
    reqs, tags = [], []

    def send(expr, tag=None):
        reqs.append({"id": len(reqs), "expr": expr})
        tags.append(tag)

    for b in benches:
        if b["setup"]:
            send(b["setup"])
        # warm-up, untimed: first touch pays page faults and any one-off compile.
        # Skipped for a `cold` benchmark: some builtins memoise so hard that no
        # amount of argument variation defeats them (Mathilda's Pi cache keeps
        # guard digits, so asking for 100003 after 100000 is a cache hit and the
        # timed reps read 1 us). For those the FIRST, cold call is the honest
        # measurement of the work, so it is the one reported.
        if not b["cold"]:
            send("kk = 0; " + b["expr"])
        # {min, max} over reps. The MAX is the caching tripwire: a builtin with an
        # internal cache (PrimePi's table, Pi's digits) answers the 2nd rep in
        # ~0 s, and reporting the min alone would make the row vacuous. `kk` is
        # bound so an expression can vary its argument per rep and defeat that.
        send("tsq = Table[First[AbsoluteTiming[%s]],{kk,1,%d}]; {Min[tsq], Max[tsq]}"
             % (b["expr"], 1 if b["cold"] else b["reps"]), ("t", b["id"]))
        if b["check"]:
            send("InputForm[%s]" % b["check"], ("v", b["id"]))
    send("Quit[]")

    payload = "".join(json.dumps(r) + "\n" for r in reqs) + \
              json.dumps({"type": "quit"}) + "\n"
    t0 = time.time()
    p = subprocess.run([MATHILDA], input=payload, capture_output=True, text=True)
    out = {}
    for ln in p.stdout.splitlines():
        ln = ln.strip()
        if not ln.startswith("{"):
            continue
        try:
            o = json.loads(ln)
        except Exception:
            continue
        if o.get("type") != "expr":
            continue
        tag = tags[o["id"]] if o["id"] < len(tags) else None
        if tag:
            out.setdefault(tag[1], {})[tag[0]] = o.get("payload", "").strip()
    sys.stderr.write("[mathilda] %.1fs wall\n" % (time.time() - t0))
    if p.returncode != 0:
        sys.stderr.write("[mathilda] exit %d\n%s\n" % (p.returncode, p.stderr[-2000:]))
    return out


# --------------------------------------------------------------------------
# Mathematica: one wolframscript file (startup is ~3 s, so batch everything).
# --------------------------------------------------------------------------
def run_wolfram(benches, tmpdir):
    lines = ["$HistoryLength = 0;"]
    for b in benches:
        setup = b["wl_setup"] or b["setup"]
        expr = b["wl"] or b["expr"]
        if setup:
            lines.append(setup)
        if not b["cold"]:
            lines.append("kk = 0; " + expr)                    # warm-up
        lines.append('tsq = Table[First[AbsoluteTiming[%s]],{kk,1,%d}]; '
                     'Print["@t|%s|", {Min[tsq], Max[tsq]}];'
                     % (expr, 1 if b["cold"] else b["reps"], b["id"]))
        if b["check"]:
            chk = b["check"] if not b["wl"] else b["check"]
            lines.append('Print["@v|%s|", InputForm[%s]];' % (b["id"], chk))
    path = os.path.join(tmpdir, "hpc_bench.wls")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    t0 = time.time()
    p = subprocess.run([WOLFRAM, "-file", path], capture_output=True, text=True)
    out = {}
    for ln in p.stdout.splitlines():
        m = re.match(r"@([tv])\|([^|]+)\|(.*)", ln.strip())
        if m:
            v = m.group(3).strip()
            # Print[InputForm[x]] emits the literal wrapper "InputForm[...]" in
            # wolframscript, so every value comparison read as a mismatch. Strip it.
            if v.startswith("InputForm[") and v.endswith("]"):
                v = v[len("InputForm["):-1]
            out.setdefault(m.group(2), {})[m.group(1)] = v
    sys.stderr.write("[wolfram] %.1fs wall\n" % (time.time() - t0))
    if p.stderr.strip():
        sys.stderr.write("[wolfram] stderr:\n" + p.stderr[-2000:] + "\n")
    return out


# --------------------------------------------------------------------------
# NumPy: one generated Python file, run under an interpreter that has numpy.
#
# Same protocol as the other two -- setup once, one untimed warm-up, then the
# minimum of `reps` timed runs, with `kk` bound per rep so a benchmark can vary
# its argument. The generated script owns a single globals dict, so state carries
# from one benchmark to the next exactly as it does in the two kernel sessions.
# --------------------------------------------------------------------------
NUMPY_PRELUDE = '''\
import math, sys, time
import numpy as np
try:
    from scipy import signal as spsig
except Exception:
    spsig = None

# NOT `sig`: the Black-Scholes benchmark binds `sig` to the volatility, and every
# benchmark shares this one globals dict, so the shorter name silently turned
# scipy.signal into 0.2 for every row that ran after it.
G = {"np": np, "math": math, "sys": sys, "spsig": spsig}


def nest(f, x, m):
    """Nest[f, x, m] -- used by several benchmarks, so it lives in the prelude."""
    for _ in range(m):
        x = f(x)
    return x


G["nest"] = nest


def _fmt(v):
    """Print a value the way the other two systems print theirs.

    numpy scalars repr as np.float64(1.5), which no cross-system value
    comparison can parse, so unwrap them first.
    """
    if isinstance(v, np.generic):
        v = v.item()
    if isinstance(v, bool):
        return "True" if v else "False"
    if isinstance(v, float):
        return repr(v)
    if isinstance(v, int):
        return str(v)
    return str(v)


def _run(bid, setup, expr, chk, reps):
    try:
        if setup:
            exec(setup, G)
        code = compile(expr, "<" + bid + ">", "eval")
        G["kk"] = 0
        eval(code, G)                                    # warm-up, untimed
        ts = []
        for k in range(1, reps + 1):
            G["kk"] = k
            t0 = time.perf_counter()
            eval(code, G)
            ts.append(time.perf_counter() - t0)
        print("@t|%s|{%r, %r}" % (bid, min(ts), max(ts)), flush=True)
        if chk:
            print("@v|%s|%s" % (bid, _fmt(eval(chk, G))), flush=True)
    except Exception as e:
        sys.stderr.write("[numpy] %s FAILED: %r\\n" % (bid, e))
'''


def run_numpy(benches, tmpdir):
    lines = [NUMPY_PRELUDE]
    for b in benches:
        if not b.get("py"):
            continue
        lines.append("_run(%s, %s, %s, %s, %d)" % (
            json.dumps(b["id"]), json.dumps(b.get("py_setup") or ""),
            json.dumps(b["py"]), json.dumps(b.get("py_check") or ""),
            1 if b["cold"] else b["reps"]))
    path = os.path.join(tmpdir, "hpc_bench_np.py")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    t0 = time.time()
    p = subprocess.run([PYTHON, path], capture_output=True, text=True)
    out = {}
    for ln in p.stdout.splitlines():
        m = re.match(r"@([tv])\|([^|]+)\|(.*)", ln.strip())
        if m:
            out.setdefault(m.group(2), {})[m.group(1)] = m.group(3).strip()
    sys.stderr.write("[numpy] %.1fs wall\n" % (time.time() - t0))
    if p.stderr.strip():
        sys.stderr.write("[numpy] stderr:\n" + p.stderr[-3000:] + "\n")
    return out


def values_agree(a, b):
    """Do two printed values agree?

    Exact match first, then a numeric comparison to the SHORTER of the two
    printings: Mathilda's printer emits 6 significant figures where Wolfram emits
    17, so -5317.09 and -5317.087608456237 are the same answer and a string
    comparison would call them different.
    """
    if a == b:
        return True
    try:
        x, y = float(a.replace("*^", "e")), float(b.replace("*^", "e"))
    except ValueError:
        return False
    if x == y:
        return True
    # Significant digits of the SHORTER printing: strip sign, decimal point and
    # leading zeros (0.123564 has six significant digits, not seven), then allow
    # half a unit in the last of them, since printing rounds.
    def sigdigits(t):
        d = t.replace("-", "").replace(".", "").lstrip("0")
        return len(d.rstrip("0")) or 1
    sig = min(sigdigits(a), sigdigits(b))
    tol = 5.0 * 10.0 ** -sig
    return abs(x - y) <= tol * max(abs(x), abs(y), 1e-300)


def as_pair(s):
    """Parse a printed {min, max} into a (min, max) float pair."""
    if not s:
        return (None, None)
    t = s.replace("*^", "e").replace("`", "").strip().strip("{}")
    parts = [x.strip() for x in t.split(",")]
    try:
        v = [float(x) for x in parts]
    except ValueError:
        return (None, None)
    if len(v) == 1:
        return (v[0], v[0])
    return (v[0], v[1])


def fmt(t):
    if t is None:
        return "—"
    if t < 1e-3:
        return "%.0f us" % (t * 1e6)
    if t < 1.0:
        return "%.2f ms" % (t * 1e3)
    return "%.3f s" % t


def ratio(other, mine):
    """How Mathilda reads against another system: '3.10x' = Mathilda is faster."""
    if not other or not mine:
        return "—"
    return ("%.2fx" % (other / mine)) if other >= mine else ("1/%.2fx" % (mine / other))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--only", default="")
    ap.add_argument("--system", default="all",
                    help="comma-separated: mathilda, wolfram, numpy (default all)")
    ap.add_argument("--json", default="")
    a = ap.parse_args()

    want = set(x.strip() for x in a.system.split(",") if x.strip())
    if "all" in want or "both" in want:
        want |= {"mathilda", "wolfram", "numpy"}

    sel = set(x.strip() for x in a.only.split(",") if x.strip())
    benches = [apply_scale(b, a.scale) for b in B if not sel or b["id"] in sel]
    tmp = os.environ.get("CLAUDE_JOB_DIR", "/tmp")
    tmp = os.path.join(tmp, "tmp") if os.path.isdir(os.path.join(tmp, "tmp")) else tmp

    md = run_mathilda(benches) if "mathilda" in want else {}
    wl = run_wolfram(benches, tmp) if "wolfram" in want else {}
    npy = run_numpy(benches, tmp) if "numpy" in want else {}

    rows, group = [], None
    print("| Benchmark | Mathilda | Mathematica 14.0 | NumPy 2.4.4 | vs WL | vs NumPy |")
    print("|---|---:|---:|---:|---:|---:|")
    for b in benches:
        i = b["id"]
        (tm, tmx) = as_pair(md.get(i, {}).get("t"))
        (tw, twx) = as_pair(wl.get(i, {}).get("t"))
        (tp, tpx) = as_pair(npy.get(i, {}).get("t"))
        vm = md.get(i, {}).get("v")
        vw = wl.get(i, {}).get("v")
        vp = npy.get(i, {}).get("v")
        if b["group"] != group:
            group = b["group"]
            print("| **%s** | | | | | |" % group)
        flag = ""
        # Every system that produced a check value is compared against Mathilda;
        # a benchmark whose three columns do not agree is not timing one thing.
        if b["check"] and vm:
            if vw and not values_agree(vm, vw):
                flag += "  ⚠ WL VALUE MISMATCH"
            if vp and not values_agree(vm, vp):
                flag += "  ⚠ NUMPY VALUE MISMATCH"
        # min far below max => an internal cache answered the later reps
        for who, lo, hi in (("mathilda", tm, tmx), ("wolfram", tw, twx),
                            ("numpy", tp, tpx)):
            if lo and hi and lo > 0 and hi / lo > 5:
                flag += f"  ⚠ {who} max/min={hi/lo:.0f}x (cache?)"
        print("| %s | %s | %s | %s | %s | %s |%s"
              % (b["name"], fmt(tm), fmt(tw), fmt(tp),
                 ratio(tw, tm), ratio(tp, tm), flag))
        rows.append(dict(id=i, name=b["name"], group=b["group"],
                         mathilda=tm, mathilda_max=tmx,
                         wolfram=tw, wolfram_max=twx,
                         numpy=tp, numpy_max=tpx,
                         vm=vm, vw=vw, vp=vp,
                         # values_agree, NOT `vm == vw`: Mathilda prints 6
                         # significant figures where Wolfram prints 17, so every
                         # inexact row failed the exit check while the table
                         # beside it -- which already used values_agree -- showed
                         # no mismatch flag. The two disagreed about the same
                         # numbers.
                         agree=(None if not (b["check"] and vm)
                                else (bool(not vw or values_agree(vm, vw))
                                      and bool(not vp or values_agree(vm, vp))))))
    if a.json:
        with open(a.json, "w") as f:
            json.dump(rows, f, indent=1)
    bad = [r for r in rows if r["agree"] is False]
    if bad:
        sys.stderr.write("\nVALUE MISMATCHES:\n")
        for r in bad:
            sys.stderr.write("  %-12s mathilda=%s  wolfram=%s  numpy=%s\n"
                             % (r["id"], r["vm"], r["vw"], r["vp"]))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
