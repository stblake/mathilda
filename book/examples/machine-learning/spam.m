# Messages become numbers: exclamation marks, capital letters, and money words per message.
features[s_] := {StringCount[s, "!"], StringCount[s, RegularExpression["[A-Z]"]]/StringLength[s], StringCount[s, RegularExpression["(?i)free|win|prize|cash|claim|urgent"]]}
spam = {"WIN a FREE prize now!!!", "URGENT: claim your cash prize", "You have won! Claim FREE cash today!", "Free entry, win big!!", "CLAIM YOUR PRIZE NOW", "Cash bonus waiting, act now!"};
ham = {"Are we still on for lunch tomorrow?", "Here are the notes from the meeting.", "Can you send me the report by Friday?", "Thanks, see you at the station at six.", "The kids loved the museum on Sunday.", "I pushed the fix, please review."};
filter = Classify[Join[Thread[Map[features, spam] -> "spam"], Thread[Map[features, ham] -> "ham"]], Method -> "NaiveBayes"];
features["Congratulations! You WIN a free cruise"]
Map[filter[features[#]] &, {"Congratulations! You WIN a free cruise", "Lunch at noon? I will bring the slides.", "Claim your free cash now!!!"}]
filter[features["Meeting moved to Friday, free parking in the lot"], "Probabilities"]
