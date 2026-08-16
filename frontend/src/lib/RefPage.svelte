<script lang="ts">
  /*
   * A symbol's reference page, rendered from the markdown that
   * `python3 site/generate.py` produces and `npm run sync:refpages` mirrors
   * into public/refpages/.
   *
   * The content is deliberately NOT authored here. It is the same page the
   * published site serves, examples included -- and those examples were
   * re-verified against the built binary at generation time, which is the whole
   * reason to render this file rather than compose something new from the
   * docstring.
   */
  import { marked } from 'marked';
  import { openUrl } from './ipc';

  /* The Markdown to render. The page is fetched and split into cells by
     openRefpage, so this component is a pure renderer -- one prose segment of a
     reference page, with the examples living in real code cells around it. */
  export let markdown: string;
  /* Called with a symbol name when the reader clicks a link to another reference page. Passed in
     rather than reached for, so this stays a renderer with no knowledge of the canvas. */
  export let onOpen: ((name: string) => void) | null = null;

  let html = '';

  /* MkDocs admonitions ("!!! success \"Status: Stable\"" followed by an indented
     body) are not markdown -- marked would render the bang line as a paragraph
     and swallow the body into a code block, because it is indented four spaces.
     Lift them out to HTML first. */
  function admonitions(md: string): string {
    const lines = md.split('\n');
    const out: string[] = [];
    for (let i = 0; i < lines.length; i++) {
      const m = /^!!!\s+(\w+)\s*(?:"([^"]*)")?\s*$/.exec(lines[i]);
      if (!m) { out.push(lines[i]); continue; }
      const [, kind, title] = m;
      const body: string[] = [];
      i++;
      while (i < lines.length && (lines[i].trim() === '' || /^\s{4}/.test(lines[i]))) {
        if (lines[i].trim() !== '') body.push(lines[i].replace(/^\s{4}/, ''));
        else if (body.length) body.push('');
        i++;
      }
      i--;                                    /* the loop's i++ re-reads this line */
      out.push(
        `<div class="adm adm-${kind}">` +
        (title ? `<div class="adm-title">${escapeHtml(title)}</div>` : '') +
        `<div class="adm-body">${marked.parse(body.join('\n').trim(), { async: false })}</div>` +
        `</div>`
      );
      /* A markdown HTML block runs until a BLANK line, so without this the
         heading immediately after the admonition gets swallowed into the raw
         block and never renders -- which ate every page's "## Description",
         since the status admonition is always followed directly by it. */
      out.push('');
    }
    return out.join('\n');
  }

  function escapeHtml(s: string): string {
    return s.replace(/[&<>"']/g, c =>
      ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c] as string));
  }

  /* Stable id for a heading, shared with the table of contents so its links
     resolve. Must match tocSlug() in refpages.ts. */
  function slug(text: string): string {
    return 'ref-' + text.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');
  }

  $: {
    /* The H1 repeats the symbol name, which is already the card title. */
    const src = markdown.replace(/^#\s+.*\n/, '');
    const raw = marked.parse(admonitions(src), { async: false }) as string;
    /* Give every heading an id. The page is split across many cells, so the
       table of contents can only jump by element id -- there is no single
       document to anchor into. */
    html = raw.replace(/<(h[23])>(.*?)<\/\1>/g,
      (_m, tag, inner) => `<${tag} id="${slug(inner.replace(/<[^>]+>/g, ''))}">${inner}</${tag}>`);
  }

  /* Links in a generated page come in three kinds, and only two were handled: a `#anchor` scrolls,
     an `http(s)` URL opens in the browser, and a RELATIVE `.md` path is another reference page --
     which fell through to preventDefault and then nothing, so every "See also" link was dead. */
  function onClick(ev: MouseEvent) {
    const a = (ev.target as HTMLElement)?.closest('a');
    const href = a?.getAttribute('href');
    if (!href) return;
    ev.preventDefault();
    if (href.startsWith('#')) {
      /* A table-of-contents link. The target heading lives in a different cell,
         so this is a scroll within the card rather than navigation. */
      document.getElementById(href.slice(1))
              ?.scrollIntoView({ behavior: 'smooth', block: 'start' });
      return;
    }
    if (/^https?:/.test(href)) { openUrl(href); return; }
    /* `../category/Name.md` (or `Name.md`) is a sibling reference page: open it here rather than
       treating it as a dead end. The symbol name is the file stem, which is exactly what
       openRefpage takes. */
    const md = href.match(/([A-Za-z$][A-Za-z0-9$]*)\.md(?:#.*)?$/);
    if (md) {
      const sym = md[1];
      if (sym !== 'index' && onOpen) { onOpen(sym); return; }
    }
  }
</script>

<!-- svelte-ignore a11y-no-static-element-interactions a11y-click-events-have-key-events -->
<div class="refpage" on:click={onClick}>
  {@html html}
</div>

<style>
  .refpage {
    padding: 0.5rem 1.1rem 1.4rem;
    font: 0.86rem/1.62 var(--sans);
    color: var(--text);
    max-width: 62rem;
    /* #app sets text-align: center for the shell's own layout, and it inherits
       all the way down here -- which centred every heading, paragraph, code
       block and bullet on the page. Prose is left-aligned. */
    text-align: left;
  }

  .ref-note { color: var(--text-dim, var(--text)); font-size: 0.84rem; padding: 0.6rem 0; }
  .ref-error { color: var(--err, #f38ba8); }
  .ref-hint { color: var(--text-dim, var(--text)); font-size: 0.78rem; margin-top: 0.35rem; line-height: 1.5; }
  .ref-detail { opacity: 0.65; }

  /* ---- Block rhythm ---------------------------------------------------- */
  .refpage :global(h2) {
    font-size: 0.95rem;
    font-weight: 650;
    color: var(--text-h);
    margin: 1.5rem 0 0.5rem;
    padding-bottom: 0.28rem;
    border-bottom: 1px solid var(--border);
    letter-spacing: 0.01em;
  }
  .refpage :global(h3) {
    font-size: 0.87rem;
    font-weight: 600;
    color: var(--text-h);
    margin: 1.1rem 0 0.4rem;
  }
  .refpage :global(h2:first-child), .refpage :global(h3:first-child) { margin-top: 0.2rem; }
  .refpage :global(p) { margin: 0.5rem 0; }
  .refpage :global(ul), .refpage :global(ol) { margin: 0.45rem 0; padding-left: 1.35rem; }
  .refpage :global(li) { margin: 0.22rem 0; }
  .refpage :global(strong) { color: var(--text-h); font-weight: 620; }
  .refpage :global(hr) { border: 0; border-top: 1px solid var(--border); margin: 1.2rem 0; }

  /* ---- Code ------------------------------------------------------------ */
  .refpage :global(code) {
    font: 0.84em/1.5 var(--mono);
    background: var(--surface-2);
    border-radius: 4px;
    padding: 0.06em 0.32em;
  }
  .refpage :global(pre) {
    background: var(--surface-2);
    border: 1px solid var(--border);
    border-radius: 7px;
    padding: 0.6rem 0.8rem;
    overflow-x: auto;               /* long example lines scroll, page does not */
    margin: 0.55rem 0;
  }
  .refpage :global(pre code) {
    background: none;
    padding: 0;
    font-size: 0.8rem;
    line-height: 1.65;
    /* The generated examples are In[n]:= / Out[n]= transcripts; keeping them
       unwrapped preserves the alignment they were written with. */
    white-space: pre;
  }

  /* ---- Tables (the Method matrix, option tables) ----------------------- */
  .refpage :global(table) {
    border-collapse: collapse;
    margin: 0.6rem 0;
    font-size: 0.8rem;
    display: block;
    overflow-x: auto;               /* wide matrices scroll in place */
    max-width: 100%;
  }
  .refpage :global(th), .refpage :global(td) {
    border: 1px solid var(--border);
    padding: 0.26rem 0.55rem;
    text-align: left;
  }
  .refpage :global(th) { background: var(--surface-2); color: var(--text-h); font-weight: 600; }

  /* ---- Links ----------------------------------------------------------- */
  .refpage :global(a) {
    color: var(--accent, #89b4fa);
    text-decoration: underline;
    text-underline-offset: 2px;
    cursor: pointer;
  }

  /* ---- Collapsed notes -------------------------------------------------- */
  .refpage :global(details) {
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 0.35rem 0.7rem;
    margin: 0.6rem 0;
    background: var(--surface-2);
  }
  .refpage :global(summary) {
    cursor: pointer;
    font-weight: 600;
    color: var(--text-h);
    font-size: 0.82rem;
    list-style-position: outside;
  }
  .refpage :global(details[open] summary) { margin-bottom: 0.3rem; }

  /* ---- Table of contents ------------------------------------------------ */
  .refpage :global(.ref-toc) { margin: 0.2rem 0 0.4rem; }
  .refpage :global(.ref-toc ul) { margin: 0.2rem 0; padding-left: 1.1rem; }
  .refpage :global(.ref-toc li) { margin: 0.1rem 0; }

  /* ---- Admonitions ----------------------------------------------------- */
  .refpage :global(.adm) {
    border: 1px solid var(--border);
    border-left-width: 3px;
    border-radius: 6px;
    padding: 0.5rem 0.75rem;
    margin: 0.7rem 0;
    background: var(--surface-2);
    font-size: 0.82rem;
  }
  .refpage :global(.adm-title) { font-weight: 620; color: var(--text-h); margin-bottom: 0.15rem; }
  .refpage :global(.adm-body > :first-child) { margin-top: 0; }
  .refpage :global(.adm-body > :last-child) { margin-bottom: 0; }
  .refpage :global(.adm-success) { border-left-color: var(--ok, #4ade80); }
  .refpage :global(.adm-warning) { border-left-color: #f9e2af; }
  .refpage :global(.adm-danger), .refpage :global(.adm-failure) { border-left-color: var(--err, #f38ba8); }
  .refpage :global(.adm-note), .refpage :global(.adm-info) { border-left-color: var(--accent, #89b4fa); }
</style>
