# 9.5 Writing and reading whole expressions: Write, Put, Get, and >>.
# Write renders each argument in input form and appends a newline -- and it
# EVALUATES the argument first. WriteString instead writes raw text verbatim.
s = OpenWrite["/tmp/mathilda_book_wr.txt"];
Write[s, 2 + 2];
Write[s, x^2 + 1];
WriteString[s, "-- a raw line --\n"];
Close[s];
ReadList["/tmp/mathilda_book_wr.txt", String]
# The >> operator (Put) writes one expression, truncating the file.
FactorInteger[40320] >> "/tmp/mathilda_book_facs.m"
Get["/tmp/mathilda_book_facs.m"]
# >>> (PutAppend) adds more without erasing what is there.
FactorInteger[12] >>> "/tmp/mathilda_book_facs.m"
FactorInteger[100] >>> "/tmp/mathilda_book_facs.m"
# Get evaluates each expression and returns the LAST; ReadList collects ALL.
Get["/tmp/mathilda_book_facs.m"]
ReadList["/tmp/mathilda_book_facs.m"]
# Put writes several expressions at once.
Put[a^2 + b, Sin[x], "/tmp/mathilda_book_exprs.m"]
ReadList["/tmp/mathilda_book_exprs.m"]
