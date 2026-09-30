# Featured: inverse kinematics. Choose the two joint angles of an arm with links 1 and 0.7
# so the hand reaches the target (1.2, 0.8).
hand = {Cos[a] + 0.7 Cos[a + b], Sin[a] + 0.7 Sin[a + b]};
sol = FindMinimum[Total[(hand - {1.2, 0.8})^2], {{a, 0.2}, {b, 0.4}}]
elbow = {Cos[a], Sin[a]} /. Last[sol]; tip = hand /. Last[sol];
fig = Graphics[{Thick, Blue, Line[{{0, 0}, elbow, tip}], Black, Disk[{0, 0}, 0.05], Disk[elbow, 0.04], Red, Disk[{1.2, 0.8}, 0.05]}, Frame -> True, AspectRatio -> Automatic, PlotRange -> {{-0.2, 1.6}, {-0.2, 1.2}}]
