#!/bin/bash
# Render every experiment write-up to PDF with pandoc.
#
#     docs/experiments/build_pdfs.sh              # all of them
#     docs/experiments/build_pdfs.sh 12-graph-and-sparse
#
# Each <NN>-<slug>/README.md becomes <NN>-<slug>/<slug>.pdf, and the index
# becomes experiments.pdf.  The PDFs are build products: they are regenerated
# from the markdown and are not the source of anything.
#
# NOTES ON THE RENDERING, each of which was a failure before it was a setting:
#
#   --pdf-engine=xelatex   the write-ups use × ≈ ≥ µ ⁷ ² and box-drawing in
#                          places; pdflatex cannot encode those and dies with
#                          an "Unicode character not set up" error.
#   -V mainfont            xelatex's default has no glyphs for × or µ, and
#                          Helvetica Neue BOLD has no → -- which only shows up
#                          in headings and bold table cells, so it is easy to
#                          miss. Arial Unicode MS covers everything the
#                          write-ups use; override with PDF_MAINFONT if you
#                          prefer something with more character and are willing
#                          to check the warnings.
#   -V monofont            code blocks hold the same characters.
#   --from=gfm+...         the tables are GitHub-flavoured pipe tables and the
#                          footnote markers (ᵃ) are literal superscripts, not
#                          pandoc footnotes.
#   -V colorlinks          the default link colouring draws boxes that survive
#                          into print.
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXP="$ROOT/docs/experiments"

command -v pandoc >/dev/null 2>&1 || { echo "pandoc not found" >&2; exit 1; }

ENGINE="${PDF_ENGINE:-xelatex}"
command -v "$ENGINE" >/dev/null 2>&1 || { echo "$ENGINE not found" >&2; exit 1; }

# A font with the mathematical and superscript glyphs the write-ups use.
MAIN="${PDF_MAINFONT:-Arial Unicode MS}"
MONO="${PDF_MONOFONT:-Menlo}"

render() {          # render <input.md> <output.pdf> <title>
  pandoc "$1" -o "$2" \
    --from=gfm+tex_math_dollars \
    --pdf-engine="$ENGINE" \
    --toc --toc-depth=2 \
    -V geometry:margin=2.2cm \
    -V mainfont="$MAIN" \
    -V monofont="$MONO" \
    -V fontsize=10pt \
    -V colorlinks=true \
    -V linkcolor=RoyalBlue \
    -V urlcolor=RoyalBlue \
    -V title="$3" \
    -V date="Mathilda — execution-speed experiments" \
    2>&1 | grep -vE "^\\[WARNING\\] (Could not fetch|Deprecated)" || true
}

if [ $# -ge 1 ]; then
  DIRS="$EXP/$1"
else
  DIRS=$(ls -d "$EXP"/[0-9][0-9]-* 2>/dev/null)
  # The index gets its own PDF, named for what it is rather than for a number.
  echo "  experiments.pdf"
  render "$EXP/README.md" "$EXP/experiments.pdf" "Execution-speed experiments"
fi

for d in $DIRS; do
  [ -f "$d/README.md" ] || continue
  slug=$(basename "$d" | sed 's/^[0-9][0-9]-//')
  num=$(basename "$d" | sed 's/-.*//')
  # The document title is the write-up's own H1, minus the leading "# ".
  title=$(head -1 "$d/README.md" | sed 's/^# *//')
  echo "  $(basename "$d")/$slug.pdf"
  render "$d/README.md" "$d/$slug.pdf" "$title"
done
