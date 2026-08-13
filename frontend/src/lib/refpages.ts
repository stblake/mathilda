/*
 * The reference-page index: which symbols have documentation, and where.
 *
 * Shared by the F1 hotkey and by RefPage.svelte. The hotkey needs a SYNCHRONOUS
 * answer -- a CodeMirror keybinding must return true/false immediately to decide
 * whether it handled the key -- so the index is fetched once at module load and
 * the lookup reads the resolved copy.
 */

let index: Record<string, string> | null = null;
let pending: Promise<Record<string, string>> | null = null;

/** Fetch (once) and cache the name -> path index. */
export function loadRefpageIndex(): Promise<Record<string, string>> {
  if (index) return Promise.resolve(index);
  if (!pending) {
    pending = fetch('/refpages/index.json')
      .then(r => {
        if (!r.ok) throw new Error(`index ${r.status}`);
        return r.json();
      })
      .then((data: Record<string, string>) => (index = data))
      .catch(err => {
        pending = null;              /* allow a later retry */
        throw err;
      });
  }
  return pending;
}

/** A symbol's reference page as Markdown. Throws if it has none.
 *
 * `no-store`: the pages are regenerated whenever the docs or the binary change,
 * and a cached copy silently serves the previous build -- which looked exactly
 * like the generator having failed to apply a change. */
export async function fetchRefpageMarkdown(name: string): Promise<string> {
  const index = await loadRefpageIndex();
  const rel = index[name];
  if (!rel) throw new Error(`no reference page for ${name}`);
  const res = await fetch(`/refpages/${rel}`, { cache: 'no-store' });
  if (!res.ok) throw new Error(`${rel} ${res.status}`);
  return res.text();
}

/** Saved figures for a page, keyed by the example's input text.
 *
 * Generated alongside the pages so a graphics example can show its plot without
 * the reader running it. Absent for most symbols, and absent for individual
 * figures too large to ship, in which case the example still runs. */
export async function fetchRefpageFigures(name: string): Promise<Record<string, object>> {
  const index = await loadRefpageIndex();
  const rel = index[name];
  if (!rel) return {};
  try {
    const res = await fetch(`/refpages/plots/${rel.replace(/\.md$/, '.json')}`,
                            { cache: 'no-store' });
    return res.ok ? await res.json() : {};
  } catch {
    return {};
  }
}

/** A reference page cut into renderable prose and runnable examples. */
export type RefSegment =
  | { kind: 'md'; text: string }
  | { kind: 'heading'; level: 2 | 3; text: string }
  | { kind: 'example'; input: string; output: string };

/** Split a page so its examples become real cells.
 *
 * The generated pages carry examples as In[n]:= / Out[n]= transcripts inside
 * ```mathematica fences. Rendered as a code block they are inert text; split
 * out, each becomes an input the reader can edit and re-run, seeded with the
 * output the generator verified against the built binary. So the page looks
 * complete before anything is run, and stays honest after.
 *
 * Only mathematica fences are split. The ```text fences hold ASCII tables and
 * algorithm prose, which are not expressions and must survive as written. */
export function splitRefpage(md: string): RefSegment[] {
  const out: RefSegment[] = [];
  const prose: string[] = [];
  const lines = md.split('\n');

  const flushProse = () => {
    const text = prose.join('\n').trim();
    /* Must contain something a reader can see. Trimming alone is not enough:
       stripping fenced examples and empty headings out of a section can leave
       stray punctuation -- a lone ';' or '.' -- and each such segment became a
       cell that rendered as nothing but still occupied a row, which is what
       produced the long blank gaps under a collapsed section. */
    if (/[A-Za-z0-9]/.test(text)) out.push({ kind: 'md', text });
    prose.length = 0;
  };

  for (let i = 0; i < lines.length; i++) {
    /* Headings become their own cells so the notebook's existing section
       folding applies to them: a section or subsection cell collapses the rows
       beneath it, which is exactly "collapsible subsections" and needs no new
       mechanism. Left inside a Markdown cell they would only ever be text. */
    const head = /^(#{2,3})\s+(.*)$/.exec(lines[i]);
    if (head) {
      flushProse();
      out.push({ kind: 'heading', level: head[1].length as 2 | 3, text: head[2].trim() });
      continue;
    }

    const fence = /^```(\w*)\s*$/.exec(lines[i]);
    if (!fence) { prose.push(lines[i]); continue; }

    /* Collect the fence body. */
    const lang = fence[1];
    const body: string[] = [];
    i++;
    while (i < lines.length && !/^```\s*$/.test(lines[i])) body.push(lines[i++]);

    if (lang !== 'mathematica' || !body.some(l => /^In\[\d+\]:=/.test(l))) {
      prose.push('```' + lang, ...body, '```');
      continue;
    }

    flushProse();
    let j = 0;
    while (j < body.length) {
      const m = /^In\[\d+\]:=\s?(.*)$/.exec(body[j]);
      if (!m) { j++; continue; }
      const input = [m[1]];
      j++;
      /* Continuation lines of a multi-line input. */
      while (j < body.length && !/^(In|Out)\[\d+\]/.test(body[j]) && body[j].trim()) {
        input.push(body[j++]);
      }
      const outLines: string[] = [];
      const om = j < body.length ? /^Out\[\d+\]=\s?(.*)$/.exec(body[j]) : null;
      if (om) {
        outLines.push(om[1]);
        j++;
        while (j < body.length && !/^(In|Out)\[\d+\]/.test(body[j]) && body[j].trim()) {
          outLines.push(body[j++]);
        }
      }
      out.push({
        kind: 'example',
        input: input.join('\n').trim(),
        output: outLines.join('\n').trim(),
      });
    }
  }
  flushProse();
  return out;
}

/** Must match slug() in RefPage.svelte. */
function tocSlug(text: string): string {
  return 'ref-' + text.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');
}

/** A contents list for the page, or null when there is too little to index.
 *
 * Built from the H2 headings, with H3s nested under them. The links are plain
 * fragment links; RefPage intercepts them and scrolls, because the page is
 * spread over many cells and there is no document to anchor into. */
export function buildToc(md: string): string | null {
  const items: { level: number; text: string }[] = [];
  let inFence = false;
  for (const line of md.split('\n')) {
    if (/^```/.test(line)) { inFence = !inFence; continue; }
    if (inFence) continue;                       /* '## ' inside a code block */
    const m = /^(#{2,3})\s+(.*)$/.exec(line);
    if (m) items.push({ level: m[1].length, text: m[2].trim() });
  }
  if (items.filter(i => i.level === 2).length < 3) return null;

  /* A disclosure rather than a plain list: the contents of a long page are
     worth having, but not worth pushing the first example below the fold. */
  const out = ['<details class="ref-toc" open>', '<summary>Contents</summary>', ''];
  let open = false;
  for (const it of items) {
    const link = `[${it.text}](#${tocSlug(it.text)})`;
    if (it.level === 2) {
      out.push(`- ${link}`);
      open = true;
    } else if (open) {
      out.push(`    - ${link}`);
    }
  }
  out.push('', '</details>');
  return out.join('\n');
}

/** Path of a symbol's page, or null if it has none / the index is not in yet. */
export function refpagePath(name: string): string | null {
  return index?.[name] ?? null;
}

/** Does this symbol have a reference page?
 *
 * Synchronous, and deliberately answers `false` while the index is still in
 * flight: the F1 handler uses this to decide whether to consume the keypress,
 * and consuming it to then show nothing is worse than letting it fall through.
 * The fetch is kicked off at import time below, so the window is a few ms at
 * startup rather than something a user is likely to hit. */
export function hasRefpage(name: string): boolean {
  return index != null && name in index;
}

/** True once the index has arrived, for callers that want to distinguish
 *  "no such symbol" from "not loaded yet". */
export function refpageIndexReady(): boolean {
  return index != null;
}

/* Warm the cache immediately so the first F1 press has an answer. A failure here
   is not fatal: hasRefpage() stays false, F1 falls through, and RefPage still
   surfaces the error if a page is opened another way. */
void loadRefpageIndex().catch(() => {});

/* ---------------------------------------------------------------------------
 * Cmd/Ctrl+click anywhere to open a symbol's reference page.
 *
 * Registered on the document, at module scope, deliberately. The first attempt
 * put this inside CodeMirror's extension array, which fails in a way that looks
 * exactly like "the binding does not work": an EditorView is built once in
 * CodeCell's onMount, so every cell that already existed kept an editor
 * constructed from the previous extensions and never saw the handler. A
 * document-level listener has no such lifecycle -- it is installed when this
 * module loads and covers every cell, old or new, plus output text.
 *
 * Capture phase, so it runs before CodeMirror's own mouse handling turns the
 * click into a selection.
 * ------------------------------------------------------------------------- */

/** The identifier under a viewport point, using the DOM rather than any editor
 *  API so it works over code cells and rendered output alike. */
function identifierAtPoint(x: number, y: number): string | null {
  const doc = document as Document & {
    caretRangeFromPoint?: (x: number, y: number) => Range | null;
    caretPositionFromPoint?: (x: number, y: number) => { offsetNode: Node; offset: number } | null;
  };

  let node: Node | null = null;
  let offset = 0;
  if (doc.caretRangeFromPoint) {                    /* WebKit, which Tauri uses */
    const r = doc.caretRangeFromPoint(x, y);
    if (r) { node = r.startContainer; offset = r.startOffset; }
  } else if (doc.caretPositionFromPoint) {          /* Firefox */
    const p = doc.caretPositionFromPoint(x, y);
    if (p) { node = p.offsetNode; offset = p.offset; }
  }
  if (!node || node.nodeType !== Node.TEXT_NODE) return null;

  return identifierAt(node.textContent ?? '', offset);
}

/** The identifier spanning `offset` in `text`, if any. Offset may sit just past
 *  the end of a name, which is where a caret lands after typing one. */
function identifierAt(text: string, offset: number): string | null {
  for (const m of text.matchAll(/[A-Za-z$][A-Za-z0-9$]*/g)) {
    const start = m.index ?? 0;
    const end = start + m[0].length;
    if (offset >= start && offset <= end) return m[0];
  }
  return null;
}

/* A module-level guard, because a hot reload re-runs this file and would
   otherwise stack a second identical listener on the document. */
declare global {
  interface Window {
    __mathildaDocsClick?: (e: MouseEvent) => void;
    __mathildaDocsKey?: (e: KeyboardEvent) => void;
  }
}

if (typeof document !== 'undefined') {
  if (window.__mathildaDocsClick) {
    document.removeEventListener('mousedown', window.__mathildaDocsClick, true);
  }
  const handler = (e: MouseEvent) => {
    if (!e.metaKey && !e.ctrlKey) return;
    const name = identifierAtPoint(e.clientX, e.clientY);
    /* Only claim the click for a symbol that actually has documentation, so
       Cmd+clicking a variable or prose behaves as it always did. */
    if (!name || !hasRefpage(name)) return;
    e.preventDefault();
    e.stopPropagation();
    (e.target as HTMLElement).dispatchEvent(new CustomEvent('mathilda-refpage',
      { detail: { name, at: { x: e.clientX, y: e.clientY } }, bubbles: true }));
  };
  window.__mathildaDocsClick = handler;
  document.addEventListener('mousedown', handler, true);

  /* The same gesture from the keyboard. This lives here rather than in
     CodeMirror's keymap for the reason the click handler does: a keymap only
     reaches editors constructed after it was added, so Cmd+I worked in new
     cells and silently did nothing in every cell that already existed. The
     caret position comes from the DOM selection, which is editor-agnostic. */
  if (window.__mathildaDocsKey) {
    document.removeEventListener('keydown', window.__mathildaDocsKey, true);
  }
  const keyHandler = (e: KeyboardEvent) => {
    const isI = (e.metaKey || e.ctrlKey) && (e.key.toLowerCase() === 'i' || e.code === 'KeyI');
    if (!(e.key === 'F1' || isI)) return;

    const sel = window.getSelection();
    if (!sel || !sel.rangeCount) return;

    let name: string | null = null;
    let at: { x: number; y: number } | null = null;

    const picked = sel.toString().trim();
    if (/^[A-Za-z$][A-Za-z0-9$]*$/.test(picked)) {
      name = picked;                                /* an explicit selection wins */
    }

    /* Locate the caret on screen and read the identifier there, reusing the
       lookup the click gesture uses. Reading sel.anchorNode directly is not
       enough: in CodeMirror the caret's anchor is often the line ELEMENT rather
       than a text node, with anchorOffset a child index, and the text-node test
       then failed silently -- which is why this key appeared dead while
       Cmd+click worked. The caret sits BETWEEN characters, so probe a couple of
       pixels either side before giving up. */
    const rect = sel.getRangeAt(0).getBoundingClientRect();
    if (rect.width || rect.height || rect.x || rect.y) {
      const y = rect.y + rect.height / 2;
      at = { x: rect.x, y };
      if (!name) {
        for (const dx of [-2, 2, -6, 6]) {
          name = identifierAtPoint(rect.x + dx, y);
          if (name && hasRefpage(name)) break;
        }
      }
    }

    /* Last resort: the anchor really is a text node. */
    if (!name && sel.anchorNode?.nodeType === Node.TEXT_NODE) {
      name = identifierAt(sel.anchorNode.textContent ?? '', sel.anchorOffset);
    }

    if (!name || !hasRefpage(name)) return;
    e.preventDefault();
    e.stopPropagation();
    const target = (sel.anchorNode?.parentElement ?? document.body) as HTMLElement;
    target.dispatchEvent(new CustomEvent('mathilda-refpage',
      { detail: { name, at }, bubbles: true }));
  };
  window.__mathildaDocsKey = keyHandler;
  document.addEventListener('keydown', keyHandler, true);
}
