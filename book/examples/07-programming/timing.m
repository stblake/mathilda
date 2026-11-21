# Ch.7 -- timing and pausing execution: the two clocks, Pause, $TimeUnit
$TimeUnit
AbsoluteTiming[Pause[0.3]]
Timing[Pause[0.3]]
AbsoluteTiming[Pause[1/4]]
SessionTime[]
TimeUsed[]
# TimeUsed already holds the start-up cost; difference it around a pause
t0 = TimeUsed[];
Pause[0.5]
TimeUsed[] - t0
# TimeConstrained bounds wall-clock time (v0.463): a pause is charged against the budget
AbsoluteTiming[TimeConstrained[Pause[5], 0.5]]
Timing[TimeConstrained[Pause[5], 0.5]]
TimeConstrained[Pause[5], 0.5, "too slow"]
