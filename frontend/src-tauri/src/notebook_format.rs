//! notebook_format.rs — the plain-text `.mathilda` notebook format.
//!
//! One notebook, no stored outputs, one stanza per cell:
//!
//! ```text
//! (* cell: code *)
//! Integrate[x^2, {x, 0, 1}]
//!
//! (* cell: text *)
//! Some **prose**.
//! ```
//!
//! The marker is a Mathilda comment, so a notebook of code cells is also a valid
//! script for `Mathilda -file`. Pure functions of strings, kept apart from the
//! Tauri commands so they can be unit tested without an app handle.
//!
//! What the format does NOT carry: outputs (ephemeral by design), the canvas
//! layout, and side-by-side cells -- a row of several cells is written as that
//! many consecutive stanzas.

use serde_json::{json, Value};

/// Cell types the notebook model knows. Anything else read from a file becomes
/// `code`, so a hand-edited or future file never produces a cell the UI cannot
/// render.
///
/// MUST match `CELL_STYLES` in `frontend/src/lib/notebook.ts` (plus `ref`, which is generated rather
/// than chosen). The fallback runs on SERIALIZE as well as on parse, so a style the front end offers
/// and this list omits is not merely unrecognised on load -- it is written to the file as `code` and
/// the reader's heading is gone. `npm run check:notebook` diffs the two lists for exactly that
/// reason.
const KNOWN_TYPES: &[&str] = &[
    "code",
    "text",
    "title",
    "subtitle",
    "chapter",
    "section",
    "subsection",
    "subsubsection",
    "ref",
];

/// Normalise a cell type read from a file or handed over by the front end.
/// `prose` is the pre-Markdown name for a text cell.
fn normalise_type(t: &str) -> &'static str {
    let t = t.trim();
    if t == "prose" {
        return "text";
    }
    KNOWN_TYPES
        .iter()
        .copied()
        .find(|k| *k == t)
        .unwrap_or("code")
}

/// If `line` is a stanza marker, the cell type it names.
///
/// Strict on purpose: the whole line must be `(* cell: <word> *)` (surrounding
/// whitespace aside). A looser prefix match would split a cell at any comment
/// that merely begins with "(* cell:".
fn marker_type(line: &str) -> Option<&str> {
    let inner = line
        .trim()
        .strip_prefix("(* cell:")?
        .strip_suffix("*)")?
        .trim();
    if inner.is_empty() || inner.contains(char::is_whitespace) {
        return None;
    }
    Some(inner)
}

/// Serialise cells (`[{type, source}, ...]`) to the stanza format.
/// Only `type` and `source` are read; anything else (outputs) is ignored.
pub fn serialize_stanzas(cells: &[Value]) -> String {
    let mut out = String::new();
    for (i, cell) in cells.iter().enumerate() {
        let cell_type = normalise_type(cell["type"].as_str().unwrap_or("code"));
        let source = cell["source"].as_str().unwrap_or("");
        if i > 0 {
            out.push('\n'); // blank line between stanzas
        }
        out.push_str("(* cell: ");
        out.push_str(cell_type);
        out.push_str(" *)\n");
        out.push_str(source);
        if !source.ends_with('\n') {
            out.push('\n');
        }
    }
    out
}

fn push_cell(cells: &mut Vec<Value>, cell_type: &str, source: &str) {
    cells.push(json!({
        "type": normalise_type(cell_type),
        "source": source.trim_end_matches('\n'),
    }));
}

/// Parse the stanza format into `[{type, source}, ...]`.
///
/// * Text before the first marker is kept as a leading code cell when it is not
///   blank, rather than silently dropped.
/// * A file with no markers at all is one code cell -- so an ordinary `.m`
///   script opens as a notebook.
/// * CRLF line endings are accepted (`str::lines` strips the `\r`).
pub fn parse_stanzas(content: &str) -> Vec<Value> {
    let mut cells: Vec<Value> = Vec::new();
    let mut current_type: Option<String> = None;
    let mut source = String::new();

    for line in content.lines() {
        if let Some(t) = marker_type(line) {
            match current_type.take() {
                Some(prev) => push_cell(&mut cells, &prev, &source),
                None if !source.trim().is_empty() => push_cell(&mut cells, "code", &source),
                None => {}
            }
            source.clear();
            current_type = Some(t.to_string());
        } else {
            source.push_str(line);
            source.push('\n');
        }
    }
    match current_type {
        Some(t) => push_cell(&mut cells, &t, &source),
        None if !source.trim().is_empty() => push_cell(&mut cells, "code", &source),
        None => {}
    }
    cells
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cell(t: &str, s: &str) -> Value {
        json!({ "type": t, "source": s })
    }

    #[test]
    fn serializes_readme_example() {
        let cells = vec![
            cell("code", "Integrate[x^2, {x, 0, 1}]"),
            cell("code", "Factor[x^4 - 1]"),
        ];
        assert_eq!(
            serialize_stanzas(&cells),
            "(* cell: code *)\nIntegrate[x^2, {x, 0, 1}]\n\n(* cell: code *)\nFactor[x^4 - 1]\n"
        );
    }

    #[test]
    fn round_trips_every_cell_type_and_multiline_source() {
        let cells = vec![
            cell("section", "A section"),
            cell("text", "Some **prose**\n\nwith a blank line"),
            cell("code", "a = 1\nb = 2\na + b"),
            cell("subsection", "Sub"),
            cell("code", ""),
            cell("ref", "Sin"),
        ];
        assert_eq!(parse_stanzas(&serialize_stanzas(&cells)), cells);
    }

    /// Every heading style survives a save and a reload. Worth its own test because the failure is
    /// silent and lossy in the same direction for all of them: normalise_type runs on serialize too,
    /// so a style missing from KNOWN_TYPES is written out as `code` and the heading is destroyed by
    /// the round trip rather than merely misread.
    #[test]
    fn round_trips_every_heading_style() {
        for style in ["title", "subtitle", "chapter", "section", "subsection", "subsubsection"] {
            let cells = vec![cell(style, "Heading text"), cell("code", "1 + 1")];
            assert_eq!(
                parse_stanzas(&serialize_stanzas(&cells)),
                cells,
                "{style} did not survive the round trip"
            );
        }
    }

    #[test]
    fn outputs_and_extra_fields_are_not_written() {
        let cells = vec![json!({ "type": "code", "source": "1+1", "output": [{"kind": "expr"}] })];
        let text = serialize_stanzas(&cells);
        assert!(!text.contains("output"));
        assert_eq!(parse_stanzas(&text), vec![cell("code", "1+1")]);
    }

    #[test]
    fn unknown_and_legacy_types_are_normalised() {
        let text = "(* cell: prose *)\nhello\n(* cell: widget *)\n1+1\n";
        assert_eq!(
            parse_stanzas(text),
            vec![cell("text", "hello"), cell("code", "1+1")]
        );
        assert_eq!(
            serialize_stanzas(&[cell("bogus", "x")]),
            "(* cell: code *)\nx\n"
        );
    }

    #[test]
    fn file_without_markers_is_one_code_cell() {
        assert_eq!(
            parse_stanzas("x = 2\nx^2\n"),
            vec![cell("code", "x = 2\nx^2")]
        );
        assert!(parse_stanzas("").is_empty());
        assert!(parse_stanzas("\n  \n").is_empty());
    }

    #[test]
    fn preamble_before_first_marker_is_kept() {
        let text = "(* a header comment *)\nx = 1\n\n(* cell: code *)\nx + 1\n";
        assert_eq!(
            parse_stanzas(text),
            vec![
                cell("code", "(* a header comment *)\nx = 1"),
                cell("code", "x + 1")
            ]
        );
    }

    #[test]
    fn crlf_is_accepted() {
        let text = "(* cell: code *)\r\n1+1\r\n\r\n(* cell: text *)\r\nhi\r\n";
        assert_eq!(
            parse_stanzas(text),
            vec![cell("code", "1+1"), cell("text", "hi")]
        );
    }

    #[test]
    fn marker_must_be_the_whole_line() {
        // A comment that merely starts like a marker is source, not a cell break.
        let text = "(* cell: code *)\nf[x_] := x (* cell: not a marker *)\n(* cell: two words *)\n";
        assert_eq!(
            parse_stanzas(text),
            vec![cell(
                "code",
                "f[x_] := x (* cell: not a marker *)\n(* cell: two words *)"
            )]
        );
        assert_eq!(marker_type("  (* cell: text *)  "), Some("text"));
        assert_eq!(marker_type("(* cell: *)"), None);
    }
}
