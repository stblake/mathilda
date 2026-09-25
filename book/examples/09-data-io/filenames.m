# 9.6 File names and existence.
# FileExistsQ asks the filesystem; we create a file so the test is reproducible.
s = OpenWrite["/tmp/mathilda_book_probe.txt"];
WriteString[s, "present\n"];
Close[s];
FileExistsQ["/tmp/mathilda_book_probe.txt"]
FileExistsQ["/tmp/mathilda_book_absent_98765"]
# The rest are pure string operations -- they never touch the disk.
FileExtension["report.tar.gz"]
FileBaseName["report.tar.gz"]
FileExtension["/home/sam/data.csv"]
FileBaseName["/home/sam/data.csv"]
FileNameSplit["/home/sam/mathilda/init.m"]
FileNameJoin[{"data", "2026", "run.dat"}]
FileNameJoin[FileNameSplit["/home/sam/mathilda/init.m"]]
