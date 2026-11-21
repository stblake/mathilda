# LeanVerify: the source check runs without Lean; the check by Lean does not.
r = LeanVerify["theorem ok : 2 + 2 = 4 := by decide"];
{r["Success"], r["ExitCode"], r["Replay"]}
LeanVerify["theorem nd : 2 + 2 = 4 := by native_decide"]["Rejected"]
LeanVerify["import Lean\ntheorem t : True := trivial"]["Rejected"]
LeanVerify["def f := 1\ntheorem t : f = 1 := rfl"]["Rejected"]
LeanVerify["theorem t : True := trivial -- a comment"]["Rejected"]
LeanVerify["theorem nd : 2 + 2 = 4 := by native_decide"]["ExitCode"]
