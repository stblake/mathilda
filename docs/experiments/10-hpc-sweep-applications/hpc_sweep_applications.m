(* =======================================================================
   Experiment 10 -- Second sweep: eight real-application pipelines
   =======================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md.
   
       Mathilda      ./Mathilda -file hpc_sweep_applications.m
       Mathematica   wolframscript -file hpc_sweep_applications.m
       Python        python3 hpc_sweep_applications.py
   
   WHAT IT MEASURES.  Eight pipelines out of real applications rather than eight
   builtins.  The hypothesis was that A COMPOSITION CAN BE SLOW IN WAYS NONE OF
   ITS PARTS ARE -- the previous sweep's conjugate-gradient row was slower than
   Mathematica while every primitive in it was at or ahead of parity, and only a
   composition finds that.
   
   It paid immediately: two of the eight did not finish on the first run, and the
   defect behind six of the seven fixes was the same one -- an operation with a
   working buffer path AND a working list path, quietly taking the second
   whenever the two representations met.
   
   WARNING: this file allocates several hundred megabytes.  Each kernel is sized
   as it was measured; run one at a time if memory is tight.
   ======================================================================= *)

SetAttributes[bench, HoldRest];

(* bench[label, expr] -- one untimed warm-up, then the MINIMUM of three timed
   runs.  The minimum, not the mean: we are measuring the cost of the work, and
   every source of noise on a loaded machine can only add. *)
bench[label_String, expr_] := Module[{ts},
  expr;
  ts = Table[First[AbsoluteTiming[expr]], {3}];
  Print[StringPadRight[label, 52], ToString[Round[1000. Min[ts], 0.001]], " ms"]
];

(* bench1 -- a single timed run, no warm-up, for kernels that take seconds.
   Reported as such wherever it is used. *)
SetAttributes[bench1, HoldRest];
bench1[label_String, expr_] :=
  Print[StringPadRight[label, 52],
        ToString[Round[1000. First[AbsoluteTiming[expr]], 0.001]], " ms  (1 run)"];

check[label_String, value_] :=
  Print[StringPadRight[label, 52], "check = ", value];

(* ---- bsmc -- Black-Scholes Monte Carlo, 10^7 paths --- *)
(* Computational finance: the whole option-pricing workload -- draw, Box- *)
(* Muller, lognormal, payoff, reduce, over 10^7 paths. *)
npth=10^7;
s0=100.;
strike=100.;
rr=0.05;
sig=0.2;
tt=1.;
bs[m_]:=Module[{u1,u2,z,st}, u1=RandomReal[{0,1},m]; u2=RandomReal[{0,1},m]; z=Sqrt[-2. Log[u1]] Cos[6.283185307179586 u2]; st=s0 Exp[(rr-0.5 sig^2) tt + sig Sqrt[tt] z]; Exp[-rr tt] Total[UnitStep[st-strike](st-strike)]/m];

(* ---- tseries -- Return series: EMA + rolling vol + drawdown, 10^6 -- *)
(* Sequential scans: an exponential moving average, a rolling volatility *)
(* and a maximum drawdown. *)
nts=10^6;
wts=250;
ret=RandomReal[{-0.01,0.01},nts];
tsr[r_,w_]:=Module[{px,ema,run,cs,vol}, px=Exp[Accumulate[r]]; ema=FoldList[0.98 #1+0.02 #2&,First[px],Rest[px]]; run=FoldList[Max,First[px],Rest[px]]; cs=Accumulate[r r]; vol=Sqrt[(Drop[cs,w]-Drop[cs,-w])/w]; Last[ema]+Max[1.-px/run]+Max[vol]];

(* ---- logreg -- Logistic regression, 200000x32, 100 GD steps -- *)
(* Machine learning: matrix-vector both ways with loop-carried state. *)
nlr=200000;
dlr=32;
Xlr=RandomReal[{-1,1},{nlr,dlr}];
wtr=RandomReal[{-1,1},dlr];
ylr=UnitStep[Xlr . wtr];
lrfit[xm_,yv_,m_]:=Module[{w,pr,g,k,nn,dd}, nn=Length[xm]; dd=Length[xm[[1]]]; w=ConstantArray[0.,dd]; k=0; While[k<m, pr=1./(1.+Exp[-(xm . w)]); g=(Transpose[xm] . (pr-yv))/nn; w=w-4. g; k=k+1]; Total[w w]];

(* ---- kmeans -- k-means, 100000 points, 8 dims, k = 16, 20 its -- *)
(* Masked reductions, an argmin and a scatter-as-dot. *)
nkm=100000;
dkm=8;
kkm=16;
Pkm=Transpose[RandomReal[{0,1},{nkm,dkm}]];
km[ptsT_,kc_,m_]:=Module[{cs,d,mnv,asg,mj,cnt,i,nn}, nn=Length[ptsT[[1]]]; cs=Table[ptsT[[All,j]],{j,kc}]; i=0; While[i<m, d=Table[Total[(ptsT-cs[[j]])^2],{j,kc}]; mnv=MapThread[Min,d]; asg=ConstantArray[0.,nn]; cs=Table[mj=UnitStep[mnv-d[[j]]](1.-asg); asg=asg+mj; cnt=Total[mj]; If[cnt>0.,(ptsT . mj)/cnt,cs[[j]]],{j,kc}]; i=i+1]; Total[MapThread[Min,Table[Total[(ptsT-cs[[j]])^2],{j,kc}]]]];

(* ---- nbody -- N-body all-pairs gravity, 1024 bodies, 10 steps -- *)
(* Astrophysics: Outer over all pairs, with megabyte intermediates. *)
nbd=1024;
st0=Table[RandomReal[{0,1},nbd],{6}];
grav[s_]:=Module[{px,py,pz,vx,vy,vz,dx,dy,dz,r2,iv,ax,ay,az,dt=0.001,eps=0.01}, px=s[[1]];py=s[[2]];pz=s[[3]]; vx=s[[4]];vy=s[[5]];vz=s[[6]]; dx=Outer[Subtract,px,px]; dy=Outer[Subtract,py,py]; dz=Outer[Subtract,pz,pz]; r2=dx dx+dy dy+dz dz+eps; iv=r2^(-1.5); ax=-Total[dx iv,{2}]; ay=-Total[dy iv,{2}]; az=-Total[dz iv,{2}]; {px+dt(vx+0.5 dt ax),py+dt(vy+0.5 dt ay),pz+dt(vz+0.5 dt az),vx+dt ax,vy+dt ay,vz+dt az}];

(* ---- heat3d -- 3D heat equation, 128^3, 50 steps ------- *)
(* CFD: a 3-D 7-point stencil over a 16 MB working set. *)
n3=128;
u3=RandomReal[{0,1},{n3,n3,n3}];
hstep[u_]:=u+0.1(RotateLeft[u,{1,0,0}]+RotateRight[u,{1,0,0}]+RotateLeft[u,{0,1,0}]+RotateRight[u,{0,1,0}]+RotateLeft[u,{0,0,1}]+RotateRight[u,{0,0,1}]-6. u);

(* ---- psd -- Welch PSD, 2^22 samples, 1024 blocks of 4096 -- *)
(* Signal processing: many medium FFTs, windowed. *)
nsg=2^22;
blk=4096;
nblk=1024;
sgn=RandomReal[{-1,1},nsg];
wnd=Table[0.5(1.-Cos[2. Pi (j-1)/(blk-1)]),{j,blk}];
psd[s_,w_,nb_,bl_]:=Module[{acc,i,seg}, acc=ConstantArray[0.,bl]; i=0; While[i<nb, seg=Take[s,{i bl+1,i bl+bl}] w; acc=acc+Abs[Fourier[seg]]^2; i=i+1]; acc];

(* ---- imgpipe -- Gaussian blur + Sobel edges, 1024^2 ----- *)
(* Vision: a separable Gaussian blur and a Sobel edge filter, rank-2 *)
(* correlation throughout. *)
imr=RandomReal[{0,1},{1024,1024}];
gk={{1.,4.,6.,4.,1.}}/16.;
gkT=Transpose[gk];
sxk={{-1.,0.,1.},{-2.,0.,2.},{-1.,0.,1.}};
syk={{-1.,-2.,-1.},{0.,0.,0.},{1.,2.,1.}};
edge[im_]:=Module[{bl,gx,gy}, bl=ListCorrelate[gkT,ListCorrelate[gk,im]]; gx=ListCorrelate[sxk,bl]; gy=ListCorrelate[syk,bl]; Sqrt[gx gx+gy gy]];

Print["Experiment 10 -- Second sweep: eight real-application pipelines"];
Print[""];

bench1["Black-Scholes Monte Carlo, 10^7 paths", bs[npth]];
check["Black-Scholes Monte Carlo, 10^7 paths", 100. 0.5 (1.+Erf[0.35/Sqrt[2.]]) - 100. Exp[-0.05] 0.5 (1.+Erf[0.15/Sqrt[2.]])];
bench1["Return series: EMA + rolling vol + drawdown, 10^6", tsr[ret,wts]];
check["Return series: EMA + rolling vol + drawdown, 10^6", tsr[Table[N[Mod[3 i, 7]]/700. - 0.004, {i,2000}], 50]];
bench1["Logistic regression, 200000x32, 100 GD steps", lrfit[Xlr,ylr,100]];
check["Logistic regression, 200000x32, 100 GD steps", lrfit[Table[N[Mod[i k, 5]]/5. - 0.4, {i,64},{k,4}], Table[N[Mod[i,2]], {i,64}], 50]];
bench1["k-means, 100000 points, 8 dims, k = 16, 20 its", km[Pkm,kkm,20]];
check["k-means, 100000 points, 8 dims, k = 16, 20 its", km[Transpose[Table[N[Mod[i^2+3 i k, 101]]/101.+0.001 N[k],{i,200},{k,4}]], 4, 10]];
bench1["N-body all-pairs gravity, 1024 bodies, 10 steps", Nest[grav,st0,10]];
check["N-body all-pairs gravity, 1024 bodies, 10 steps", Total[Nest[grav, Table[N[Mod[7 i + 3 k, 11]]/11., {i,6},{k,8}], 3], 2]];
bench1["3D heat equation, 128^3, 50 steps", Nest[hstep,u3,50]];
check["3D heat equation, 128^3, 50 steps", Total[Nest[hstep, Table[N[Mod[i+2j+3k,5]],{i,8},{j,8},{k,8}], 5], 3]];
bench1["Welch PSD, 2^22 samples, 1024 blocks of 4096", psd[sgn,wnd,nblk,blk]];
check["Welch PSD, 2^22 samples, 1024 blocks of 4096", Total[psd[Table[N[Sin[0.1 j]],{j,4096}], Table[0.5(1.-Cos[2. Pi (j-1)/255.]),{j,256}], 16, 256]]];
bench1["Gaussian blur + Sobel edges, 1024^2", edge[imr]];
check["Gaussian blur + Sobel edges, 1024^2", Total[edge[Table[N[Mod[i j, 17]],{i,64},{j,64}]], 2]];
