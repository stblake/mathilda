# Ch.7 -- holding evaluation: Hold, HoldForm, Evaluate, ReleaseHold
Hold[1 + 1]
HoldForm[1 + 1]
Hold[1 + 1, 2 + 2]
Hold[Evaluate[1 + 1], 2 + 2]
ReleaseHold[Hold[1 + 1]]
Attributes[Hold]
Attributes[SetDelayed]
