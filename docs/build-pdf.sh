#!/bin/sh
# Render a repository Markdown document to PDF via pandoc + xelatex.
#
#   docs/build-pdf.sh docs/compile_example/COMPILE_EXAMPLE.md [out.pdf]
#
# Run from the repository root: image paths in the documents are repo-relative.
#
# Two things need handling that a plain `pandoc -o out.pdf` gets wrong:
#
#   * <details>/<summary> is raw HTML.  The LaTeX writer DROPS raw HTML without
#     a word of complaint, so a collapsed block silently vanishes from the PDF —
#     the reader has no way to tell a section is missing.  It is rewritten into
#     an ordinary sub-heading first.
#   * Unlabelled fenced blocks become plain `verbatim`, which neither wraps nor
#     scales, so the wider bytecode listings run off the page.  The preamble
#     below shrinks verbatim and turns on breaking for highlighted blocks.
set -e

SRC=${1:?usage: docs/build-pdf.sh <input.md> [output.pdf]}
# The PDF lands beside its source unless told otherwise.  Spelled as an if
# rather than `[ test ] && OUT=...`: under `set -e` a false test there is a
# failing final command, which would abort the script for every input outside
# the current directory.
OUT=$2
if [ -z "$OUT" ]; then
    SRCDIR=$(dirname "$SRC")
    SRCBASE=$(basename "$SRC" .md)
    if [ "$SRCDIR" = "." ]; then OUT="$SRCBASE.pdf"; else OUT="$SRCDIR/$SRCBASE.pdf"; fi
fi

TMP=$(mktemp -t mathilda-pdf).md
trap 'rm -f "$TMP"' EXIT

# <summary>Foo</summary> -> a heading, so the body survives into the PDF.
# Level 3, not 4: level 4 maps to \paragraph, a RUN-IN heading, which LaTeX
# defers onto the next ordinary paragraph when a display block follows it — a
# summary followed by a code fence then prints *after* its own code.
sed -e 's|^[[:space:]]*<details>[[:space:]]*$||' \
    -e 's|^[[:space:]]*</details>[[:space:]]*$||' \
    -e 's|^[[:space:]]*<summary>\(.*\)</summary>[[:space:]]*$|### \1|' \
    "$SRC" > "$TMP"

pandoc "$TMP" -o "$OUT" \
  --from=markdown+tex_math_dollars+pipe_tables+backtick_code_blocks+fenced_code_attributes \
  --pdf-engine=xelatex \
  --resource-path=.:"$(dirname "$SRC")" \
  --highlight-style=tango \
  -V geometry:a4paper -V geometry:margin=2.2cm \
  -V fontsize=10pt \
  -V colorlinks=true -V linkcolor=RoyalBlue -V urlcolor=RoyalBlue \
  -V monofont="Menlo" \
  -H /dev/stdin <<'PREAMBLE'
\usepackage{microtype}
\usepackage{graphicx}
% Pandoc >= 3 passes the image's alt text through as \includegraphics[alt={...}],
% a key only graphicx 2022+ knows.  Accept and ignore it on older TeX Live rather
% than making authors drop the alt text (which GitHub and screen readers use).
\makeatletter
\@ifundefined{KV@Gin@alt}{\define@key{Gin}{alt}{}}{}
\makeatother
% Images: never wider than the text block.
\setkeys{Gin}{width=\linewidth,keepaspectratio}
% Code: \footnotesize is what keeps the widest bytecode listing (87 columns) on
% the page.  fvextra can additionally BREAK an over-long line instead of letting
% it run into the margin, but it is not in every TeX Live install (not in the
% 2021 basic scheme), so it is used only when present rather than required.
\fvset{fontsize=\footnotesize}
\IfFileExists{fvextra.sty}{%
  \usepackage{fvextra}%
  \DefineVerbatimEnvironment{Highlighting}{Verbatim}%
    {breaklines,breakanywhere,commandchars=\\\{\},fontsize=\footnotesize}%
}{}
% Unlabelled fenced blocks land in plain verbatim, which \fvset does not reach.
\makeatletter
\renewcommand{\verbatim@font}{\ttfamily\footnotesize}
\makeatother
% Tables in this document are wide and numeric; give them room to breathe.
\renewcommand{\arraystretch}{1.15}
% Figures stay where they were written.  In a document that says "the fourth
% panel is nearly blank", a float that drifts to another page is wrong, and
% pinning them also stops the float from forcing an overfull page.
\usepackage{float}
\floatplacement{figure}{H}
% Let a page finish short rather than be stretched to the bottom; with pinned
% figures and large code blocks, exact filling is what overflows a page.
\raggedbottom
PREAMBLE

echo "wrote $OUT ($(wc -c < "$OUT" | tr -d ' ') bytes)"
