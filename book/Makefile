# Makefile for the Mathilda book.
#
#   make pdf          regenerate links + transcripts, then build TheMathildaBook.pdf
#   make links        regenerate generated/builtinlinks.tex from builtins.json
#   make examples     run every book/examples/**/*.m through ./Mathilda
#   make check-links  fail if any \B{} names a builtin with no reference page
#   make clean        remove LaTeX build artifacts
#   make distclean    also remove the generated/ transcripts and link table
#
# The PDF build requires the Mathilda binary (../Mathilda); `make examples`
# builds it if necessary. It also requires a TeX toolchain (latexmk + pdflatex,
# biber, makeindex).

BOOK      := TheMathildaBook
PYTHON    := python3
TOOLS     := tools
GEN       := generated
MATHILDA  := ../Mathilda
LATEXMK   := latexmk
EXAMPLES  := $(shell find examples -name '*.m' 2>/dev/null)
FIGURES   := $(shell find figures -name '*.m' 2>/dev/null)

.PHONY: pdf all links version examples usage figures check-links clean distclean help

all: pdf

# --- the Mathilda binary the transcripts are generated from ---------------
$(MATHILDA):
	$(MAKE) -C .. -j

# --- generated inputs ------------------------------------------------------
links:
	$(PYTHON) $(TOOLS)/gen_links.py

version: | $(MATHILDA)
	$(PYTHON) $(TOOLS)/gen_version.py

examples: | $(MATHILDA)
	$(PYTHON) $(TOOLS)/build_examples.py

# Re-run generators whenever their inputs change (used as PDF prerequisites).
$(GEN)/builtinlinks.tex: $(TOOLS)/gen_links.py ../site/docs/assets/builtins.json
	$(PYTHON) $(TOOLS)/gen_links.py

# Re-derive the version whenever the tool changes or the binary is rebuilt (a new
# binary means a possibly new $Version, which the title/copyright pages must show).
$(GEN)/version.tex: $(TOOLS)/gen_version.py $(MATHILDA)
	$(PYTHON) $(TOOLS)/gen_version.py

$(GEN)/.examples.stamp: $(EXAMPLES) $(TOOLS)/build_examples.py ../site/verify_tutorial.py | $(MATHILDA)
	$(PYTHON) $(TOOLS)/build_examples.py
	@mkdir -p $(GEN) && touch $@

# --- docstring usage cards (\usagebox) ------------------------------------
CHAPTERS := $(shell find chapters frontmatter -name '*.tex' 2>/dev/null)
usage: | $(MATHILDA)
	$(PYTHON) $(TOOLS)/gen_usage.py

$(GEN)/.usage.stamp: $(CHAPTERS) $(TOOLS)/gen_usage.py | $(MATHILDA)
	$(PYTHON) $(TOOLS)/gen_usage.py
	@mkdir -p $(GEN) && touch $@

# --- figures (\bookfigure): plots Export-ed to PDF by the real binary ------
figures: | $(MATHILDA)
	$(PYTHON) $(TOOLS)/gen_figures.py

$(GEN)/.figures.stamp: $(FIGURES) $(TOOLS)/gen_figures.py ../site/verify_tutorial.py | $(MATHILDA)
	$(PYTHON) $(TOOLS)/gen_figures.py
	@mkdir -p $(GEN) && touch $@

# --- the book -------------------------------------------------------------
pdf: $(GEN)/builtinlinks.tex $(GEN)/version.tex $(GEN)/.examples.stamp $(GEN)/.usage.stamp $(GEN)/.figures.stamp
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error $(BOOK).tex

# --- checks ---------------------------------------------------------------
check-links:
	$(PYTHON) $(TOOLS)/check_links.py

# --- housekeeping ---------------------------------------------------------
clean:
	$(LATEXMK) -C $(BOOK).tex >/dev/null 2>&1 || true
	rm -f $(BOOK).bbl $(BOOK).run.xml $(BOOK).idx $(BOOK).ind $(BOOK).ilg *.aux chapters/*.aux frontmatter/*.aux

distclean: clean
	rm -rf $(GEN)

help:
	@sed -n '1,20p' Makefile
