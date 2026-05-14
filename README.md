# cryptography-paper-study-guides

Personal study guides for cryptography papers, written in markdown. Each guide distills a paper into a summary, the prerequisites needed to read it, its key results, and a section-by-section walkthrough.

## Layout

```
papers/
  _template/        # copy this folder when starting a new guide
  <slug>/           # one folder per paper, README.md is the guide
code/               # Rust examples that accompany select papers
```

- **Guides** live at `papers/<slug>/README.md`. The slug convention is `firstauthor-shortname-year`, e.g. `gentry-fhe-2009`.
- **Template** at `papers/_template/README.md` — copy the folder, rename, and fill in.
- **Code** lives under `code/` as a single Rust crate. Examples are added alongside the papers that motivate them; a guide's `## Code` section links to the relevant files.

## Starting a new guide

```
cp -r papers/_template papers/<slug>
$EDITOR papers/<slug>/README.md
```

## Running the code

```
cd code
cargo run
```

No examples have been added yet — the crate is a placeholder until the first paper-driven example lands.
