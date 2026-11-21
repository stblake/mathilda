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
# Working directory: move to /tmp so relative names are predictable.
FileNameTake[SetDirectory["/tmp"]]
FileNameTake[Directory[]]
s = OpenWrite["mathilda_fs_a.dat"]; WriteString[s, "1 2 3\n4 5 6\n"]; Close[s];
s = OpenWrite["mathilda_fs_b.dat"]; WriteString[s, "7 8 9\n"]; Close[s];
FileNames["mathilda_fs_*.dat"]
FileType["mathilda_fs_a.dat"]
FileType[Directory[]]
FileType["mathilda_fs_none.dat"]
DirectoryQ[ParentDirectory[]]
FileNameTake[AbsoluteFileName["mathilda_fs_a.dat"]]
ReadString["mathilda_fs_a.dat"]
ReadString["mathilda_fs_a.dat", "\n"]
FileNameTake["/home/sam/mathilda/init.m", -2]
FileNameTake["/home/sam/mathilda/init.m", {2}]
DeleteFile[FileNames["mathilda_fs_*.dat"]]
FileNames["mathilda_fs_*.dat"]
ResetDirectory[];
# A string can be read as if it were a file.
st = StringToStream["10 20 30 40"];
Skip[st, Number, 2]
Read[st, Number]
ReadString[st]
Close[st];
# Environment variables: an unset one is $Failed.
Environment["MATHILDA_BOOK_UNSET_VAR"]
StringQ[Environment["PATH"]]
