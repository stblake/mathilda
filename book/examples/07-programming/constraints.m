# Ch.7 -- MemoryConstrained and TimeConstrained share constraint frames (v0.6623)
MemoryConstrained[Range[10^7], 10^6]
MemoryConstrained[Total[Range[10^6]], 10^6, "too big"]
MemoryConstrained[RandomReal[1, 10^10], 10^3, "refused"]
MemoryConstrained[Do[Range[10^5], {100}], 10^6]
CheckAbort[MemoryConstrained[Abort[], 10^6, "budget"], "plain abort"]
TimeConstrained[CheckAbort[Pause[2], cleaned = True], 0.2]
cleaned
TimeConstrained[CheckAbort[Pause[2], "kept", PropagateAborts -> False], 0.2]
TimeConstrained[MemoryConstrained[Range[10^7], 10^6, "memory"], 5, "time"]
