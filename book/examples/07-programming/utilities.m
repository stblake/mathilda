# Ch.7 -- Echo, EchoTiming, and the parallel heads with serial semantics (v0.562)
Echo[2 + 2] + 1
Echo[{1, 2, 3}, "total: ", Total]
EchoTiming[Total[Range[10^6]]]
ParallelMap[#^2 &, {1, 2, 3}]
ParallelTable[i!, {i, 5}]
$KernelCount
