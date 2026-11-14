# AGENTS.md

Personal physics notebook: notes, derivations and simulations, for the owner's own understanding. Public repo, but not a product — keep READMEs factual and brief, no promotional language, badges or marketing.

## Structure

- `notes/` — `YYYY-MM-DD-short-topic.{md,html}` or a folder of that name with a `README.md` inside. Use the date the note is created.
- `experiments/<kebab-case-name>/` — self-contained; always include a `README.md` (what, how to run, findings/open questions).

## Conventions

- Pick the lightest tool that fits: markdown for derivations, a single self-contained HTML/JS file for interactive visualizations, Python for numerics/plots, C++ only when performance matters.
- HTML pages: one file, no build step; load libraries from a CDN with pinned versions if needed.
- Python: each experiment manages its own dependencies (e.g. a `pyproject.toml` usable with `uv`, or a `requirements.txt`). No repo-wide environment.
- C++: `CMakeLists.txt` per experiment, build into `build/` (gitignored).
- Math in markdown uses `$…$` / `$$…$$` (GitHub rendering). Use SI units unless stated otherwise; state assumptions and approximations explicitly.
- Correctness matters more than polish: when deriving, show steps; when simulating, note the numerical method and its limitations, and sanity-check against known limits or analytic results where possible.
- Don't commit generated output (large data, rendered videos, build artifacts) unless it is small and part of the note.
