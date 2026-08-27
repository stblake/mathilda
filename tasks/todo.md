# Task: Book additions — Gröbner elimination + Diophantine subsection

## Requirements
1. §4.2.7 (Gröbner bases) — add an example of **variable elimination**.
2. After §4.2.10 (Complete solutions and the reals) — new subsection on
   **solving Diophantine equations**: classical examples + algorithms discussed.
3. (follow-up) Add a **reference page for FindInstance** (regenerate site docs).
4. (follow-up) **Fully re-compile** the book incl. all examples and graphics.

## Plan
- [x] Explore book structure; identify 4.2.7 = "Gröbner bases", 4.2.10 = "Complete
      solutions and the reals" in chapters/math/algebra.tex.
- [x] Research Diophantine engine (`Solve[...,Integers]`) capabilities + algorithms.
- [x] Verify all candidate example outputs against the real binary via the NDJSON pipe.
- [x] Gröbner: appended elimination lines to examples/algebra/groebner.m (pairs 5,6:
      implicitize x=t^2,y=t^3 → x^3-y^2). Prose + Elimination Theorem callout added.
- [x] Diophantine: new examples/algebra/diophantine.m (13 codepairs) +
      diophantine-pythagorean.m + diophantine-pell.m transcripts. New subsection
      §4.2.11 ssec:alg-diophantine with 4 callouts (Hilbert 10, Pell CF, engines, Pitfall).
- [x] Reconciled forward-ref (was "§number theory") + extended section-closing paragraph.
- [x] Index concept entries added (Diophantine equation, Elimination Theorem, elimination
      ideal, Hilbert's tenth problem, Legendre's theorem, continued fraction!Pell's eqn).
- [x] Ran site/generate.py → FindInstance reference page + builtins.json entry (also
      picked up Exists/Equivalent/Implies/Xor). \B{FindInstance} now linkable.
- [x] Book recompile: links (862), examples (all verified), figures (10 PDF vector),
      pdf (85 pages, clean log, no undefined refs/citations). check-links: 345 OK.
- [x] Changelog note added under docs/spec/changelog/2026-08-24.md.

## Review / result
- Book PDF: 85 pages, TheMathildaBook.pdf. §4.2.7 gained a variable-elimination
  worked example; §4.2.11 "Diophantine equations" added after §4.2.10. Every Out[]
  produced by the real binary at build time. FindInstance reference page created.
- Enabled graphics: raylib was installed via Homebrew but the MacPorts pkg-config on
  PATH missed it; the Makefile's own PKG_CONFIG wrapper finds it, so a clean rebuild
  produced a graphics-enabled binary. 10 PDF (vector) figures regenerate headlessly.
- LIMITATION (pre-existing, environmental, chapter 3 only — NOT the algebra work):
  6 PNG raster figures (Plot3D/ComplexPlot/DensityPlot/ArrayPlot/…) crash on Export in
  a headless session — raylib PNG needs a GUI/GL context this process can't obtain.
  They degrade to graceful placeholders (\plotcell \IfFileExists). Baked in normally by
  running `make figures && make pdf` from a GUI Terminal (attached to the window server).

## Verified example outputs (from real binary)
- GroebnerBasis[{t^2-x,t^3-y},{t,x,y}] -> {x^3-y^2, -x^2+t y, t x-y, t^2-x}
- GroebnerBasis[{t^2-x,t^3-y},{x,y},{t}] -> {x^3-y^2}   (elimination ideal)
- 5x+3y==1 -> {{x->-1+3C[1], y->2-5C[1]}}; 6x+9y==5 -> {}
- {x+2y+3z==10,x-y+z==2} Integers -> {{x->18+5C[1],y->8+2C[1],z->-8-3C[1]}}
- x^2-61y^2==1 && x>0&&y>0&&x<10^10 -> {{1766319049, 226153980}}
- x^2+y^2==3z^2 -> {{0,0,0}} (Legendre: only trivial)
- x^3+y^3==1729 && 0<x<=y -> {{1,12},{9,10}} (taxicab)
- x^3+y^3+z^3==4 -> {} (mod-9 proof, no box)
- y^2==x^3-2 -> {{3,-5},{3,5}} (Mordell/Fermat)
- 2^n-7==x^2 && n>0&&x>0 -> {(3,1),(4,3),(5,5),(7,11),(15,181)} (Ramanujan-Nagell)
- n!+1==m^2 && 0<n<100 -> n=4,5,7 (Brocard)
- x^a-y^b==1 && box -> {3,2,2,3} (Catalan/Mihailescu)
- x^3+y^3==z^3 && positive -> {} (FLT/Wiles)
- x^3-2y^3==1 -> {{-1,-1},{1,0}} (Thue/Baker)
- x^2-2y^2==1 (no constraints) -> unevaluated (declines, soundness)
- x^2+y^2==z^2 && 0<x<y<z<30 -> 10 concrete triples
