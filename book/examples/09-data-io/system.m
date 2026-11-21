# System variables. Their values depend on the machine, so we test their shape.
StringQ[$System]
StringMatchQ[$SystemID, RegularExpression["[A-Za-z0-9-]+"]]
IntegerQ[$ProcessorCount] && $ProcessorCount >= 1
DirectoryQ[$InstallationDirectory]
{StringQ[$BaseDirectory], StringQ[$PreferencesDirectory]}
