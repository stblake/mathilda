# Mathilda documentation site

The source for the Mathilda documentation website, published at
<https://stblake.github.io/mathilda/>. Built with
[MkDocs Material](https://squidfunk.github.io/mkdocs-material/).

## Structure

```
site/
  mkdocs.yml          MkDocs config (theme, plugins, markdown extensions)
  requirements.txt    Python build deps (installed by CI and locally)
  generate.py         Regenerates the per-builtin reference pages
  overlays/           Hand-curated content merged onto generated pages
  docs/               The MkDocs docs_dir — what gets published
    index.md            Home page
    documentation/      Generated: one page per built-in, grouped by category
    tutorials/          Hand-written, verified tutorials
    assets/             CSS + builtins.json index
```

## How the reference pages are generated

`generate.py` builds one Markdown page per public built-in function. It:

1. Discovers builtins by parsing `symtab_set_docstring("Name", …)` calls in
   `src/**/*.c`.
2. Reads each function's attributes by driving the compiled `./Mathilda` binary.
3. Maps functions to categories from `docs/spec/builtins/*.md`.
4. Mines `In[]:=/Out[]=` examples from the spec and **re-verifies every one** by
   feeding the inputs back through `./Mathilda` and capturing the real output.
5. Derives a Stable / Partial / Experimental status.
6. Merges any hand-curated overlay from `overlays/<Name>.md`.

Because examples are verified against the binary, regeneration requires a built
`./Mathilda` (`make` at the repo root). The generated Markdown is **committed**,
so CI only needs MkDocs — no C toolchain.

### Overlay format

Drop `overlays/<Name>.md` to enrich a function's page. Front matter overrides the
auto-derived status and prepends references; the body is appended as a
"Notes & additional examples" section. **Overlay examples are not auto-verified —
verify them by hand against the binary before committing.**

```markdown
---
status: Stable
references:
  - "Author, *Title* (Publisher, year)."
---
### Worked examples

​```mathematica
In[1]:= ...
Out[1]= ...
​```
```

## Local workflow

From the repo root:

```bash
make            # build ./Mathilda (needed by the generator)
make docs       # regenerate site/docs/documentation/** from docstrings + specs

python3 -m venv site/.venv
site/.venv/bin/pip install -r site/requirements.txt
make docs-serve # live preview at http://127.0.0.1:8000
make docs-build # strict production build into site/site/
```

## Deployment

`.github/workflows/docs.yml` builds and deploys to GitHub Pages on every push to
`main` that touches `site/`. The repository's **Pages source must be set to
"GitHub Actions"** (Settings → Pages) for the workflow to publish.
