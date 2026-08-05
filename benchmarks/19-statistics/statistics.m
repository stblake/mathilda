(* Experiment 19 -- Descriptive statistics over machine arrays.
   MEASURES src/stats/ (14 files, one builtin each) through the packed path.
   Sizes are 10^6-10^7, well above PACK_MIN_ELEMENTS, so this measures buffers.
   Checks use lcgInts from data.m -- an identical integer sequence in all three
   systems, so a Mean or Median check is exactly comparable. *)

Get["../harness.m"];
Get["../data.m"];

require[{"Mean", "Median", "Variance", "StandardDeviation", "Quartiles",
         "Skewness", "Kurtosis", "RootMeanSquare", "MovingAverage",
         "MovingMedian", "ExponentialMovingAverage", "Total", "Sort",
         "Quantile", "CentralMoment"}];

n = 2000000;
v = rand01[{n}];
det = lcgInts[10000, 1000];

bench["Mean over 2x10^6", Mean[v];];
check["Mean over 2x10^6", Round[10^6 N[Mean[det]]]];

bench["Variance over 2x10^6", Variance[v];];
check["Variance over 2x10^6", Round[10^3 N[Variance[det]]]];

bench["StandardDeviation over 2x10^6", StandardDeviation[v];];
check["StandardDeviation over 2x10^6", Round[10^4 N[StandardDeviation[det]]]];

bench["Median over 2x10^6", Median[v];];
check["Median over 2x10^6", Round[10^4 N[Median[det]]]];

bench["Quartiles over 2x10^6", Quartiles[v];];
check["Quartiles over 2x10^6", Round[10^4 N[Total[Quartiles[det]]]]];

benchIf["Skewness over 2x10^6", "Skewness", Skewness[v];];
checkIf["Skewness over 2x10^6", "Skewness", Round[10^4 N[Skewness[det]]]];

benchIf["Kurtosis over 2x10^6", "Kurtosis", Kurtosis[v];];
checkIf["Kurtosis over 2x10^6", "Kurtosis", Round[10^4 N[Kurtosis[det]]]];

bench["RootMeanSquare over 2x10^6", RootMeanSquare[v];];
check["RootMeanSquare over 2x10^6", Round[10^4 N[RootMeanSquare[det]]]];

bench["MovingAverage window 100", MovingAverage[v, 100];];
check["MovingAverage window 100",
  Round[10^4 N[Total[MovingAverage[det, 100]]]]];
