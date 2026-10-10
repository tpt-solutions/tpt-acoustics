//! Checks that ticked items in `todo.md` refer to code that exists.
//!
//! This prevents the drift that left Moore loudness ticked with no
//! implementation. A ticked entry may name a Rust path in backticks, such as
//! `tpt_acou_ray::rir::hybrid_rir` (the last segment must be defined as a
//! `fn`, `struct`, `enum`, `const`, `trait`, `type` or `mod` somewhere in
//! `crates/*/src`), or a repository path such as
//! `crates/tpt-acou-ray/src/rir.rs` (the file must exist).

use std::fs;
use std::path::{Path, PathBuf};

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join("..")
}

/// Ticked entries: the `- [x]` line and its indented continuation lines.
fn ticked_entries(todo: &str) -> Vec<String> {
    let mut entries = Vec::new();
    let mut current: Option<String> = None;
    for line in todo.lines() {
        if line.starts_with("- [x]") {
            if let Some(done) = current.take() {
                entries.push(done);
            }
            current = Some(line.to_string());
        } else if line.starts_with("      ") && current.is_some() {
            current.as_mut().unwrap().push(' ');
            current.as_mut().unwrap().push_str(line.trim());
        } else {
            if let Some(done) = current.take() {
                entries.push(done);
            }
        }
    }
    if let Some(done) = current {
        entries.push(done);
    }
    entries
}

/// Backticked spans of an entry.
fn backticked(entry: &str) -> Vec<String> {
    entry
        .split('`')
        .enumerate()
        .filter(|(i, _)| i % 2 == 1)
        .map(|(_, s)| s.to_string())
        .collect()
}

fn collect_sources(dir: &Path, out: &mut Vec<PathBuf>) {
    if let Ok(entries) = fs::read_dir(dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                collect_sources(&path, out);
            } else if path.extension().is_some_and(|e| e == "rs") {
                out.push(path);
            }
        }
    }
}

fn defines(source: &str, name: &str) -> bool {
    const KEYWORDS: [&str; 7] = ["fn", "struct", "enum", "const", "trait", "type", "mod"];
    source.lines().any(|line| {
        let line = line.trim_start();
        // Strip visibility, then look for `<keyword> <name>` with a word boundary.
        let line = line
            .strip_prefix("pub(crate) ")
            .or_else(|| line.strip_prefix("pub "))
            .unwrap_or(line);
        KEYWORDS.iter().any(|kw| {
            line.strip_prefix(kw)
                .and_then(|rest| rest.strip_prefix(' '))
                .is_some_and(|rest| {
                    let ident: String = rest
                        .chars()
                        .take_while(|c| c.is_alphanumeric() || *c == '_')
                        .collect();
                    ident == name
                })
        })
    })
}

#[test]
fn ticked_todo_items_name_existing_code() {
    let root = repo_root();
    let todo = fs::read_to_string(root.join("todo.md")).expect("todo.md at the repository root");
    let mut sources = Vec::new();
    collect_sources(&root.join("crates"), &mut sources);
    let texts: Vec<String> = sources
        .iter()
        .map(|p| fs::read_to_string(p).unwrap_or_default())
        .collect();

    let mut missing: Vec<String> = Vec::new();
    for entry in ticked_entries(&todo) {
        for span in backticked(&entry) {
            if span.starts_with("crates/") {
                if !root.join(&span).exists() {
                    missing.push(format!("file {span}"));
                }
            } else if span.starts_with("tpt_acou_") && span.contains("::") {
                let name = span.rsplit("::").next().unwrap_or("");
                // Ignore method-like tails with parentheses or generics.
                let name = name.split(['(', '<']).next().unwrap_or("");
                if !name.is_empty() && !texts.iter().any(|t| defines(t, name)) {
                    missing.push(format!("symbol {span}"));
                }
            }
        }
    }
    assert!(
        missing.is_empty(),
        "ticked todo items reference code that does not exist:\n  {}",
        missing.join("\n  ")
    );
}

#[test]
fn defines_recognises_declarations() {
    assert!(defines("pub fn hybrid_rir(x: u8) {}", "hybrid_rir"));
    assert!(defines("    pub(crate) struct Tally;", "Tally"));
    assert!(!defines("let hybrid_rir_count = 3;", "hybrid_rir"));
    assert!(!defines("pub fn hybrid_rir_extra() {}", "hybrid_rir"));
}
