# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Vulgaris is a Hugo static site — a documentation-first portfolio (CV, projects, blog) for Jorge Alfredo Suazo Victoria, deployed to GitHub Pages at `https://3lconejo.github.io/Vulgaris/`.

## Commands

There is no package.json/build tooling beyond Hugo itself.

- `hugo server -D` — run the dev server with drafts included (default: http://localhost:1313/Vulgaris/)
- `hugo` — build the production site into `public/`
- `hugo --gc --minify --baseURL "<pages-url>/"` — the exact build invocation used in CI (see `.github/workflows/hugo.yml`), pinned to Hugo `0.140.2` (extended)
- `hugo new content/projects/<slug>.md` / `hugo new content/blog/<slug>.md` — scaffold new content from the archetypes in `archetypes/`

There is no test suite or linter configured. `public/` is committed (it's the built output, not gitignored) — after editing content or layouts, regenerate it with `hugo` if the built site needs to reflect changes.

## Deployment

`.github/workflows/hugo.yml` builds and deploys to GitHub Pages automatically on every push to `main`. No manual deploy step is needed.

## Architecture

Standard Hugo layout with no theme — everything is custom under `layouts/`:

- `layouts/_default/baseof.html` is the single base template all pages inherit from. It assembles `partials/head.html`, `partials/header.html`, the page's `main` block, and `partials/footer.html`, and loads `static/js/main.js`.
- `layouts/partials/` holds shared fragments: `head.html` (meta/OG tags, font/CSS links), `header.html` (nav, driven by `[menu.main]` in `hugo.toml`), `footer.html`, `card.html` (renders a single content-page preview card, used by list views), and `starfield.html` (decorative home-page background, homepage-only).
- Section-specific templates (`layouts/projects/list.html`, `layouts/blog/list.html`, `layouts/cv/single.html`, `layouts/taxonomy/`) override `_default/list.html` / `_default/single.html` for their respective content sections.
- `layouts/index.html` is the homepage template.

Content lives in `content/{projects,blog}/*.md` plus `content/cv.md`, each with frontmatter matching its archetype in `archetypes/` (`projects.md` and `blog.md` define the expected fields — e.g. projects have `status`, `summary`, `tags`, `github`, `stack`). New content should be created via `hugo new` against these archetypes rather than hand-rolled, to keep frontmatter consistent with what `card.html` and the list templates expect.

CSS/JS are plain static assets, not run through Hugo Pipes/asset processing: `static/css/{style,syntax,starfield}.css` and `static/js/main.js` are served as-is. `dev/scss/starfield.scss` is the Sass source for `static/css/starfield.css` but there is no wired-up compile step in this repo — changes to the `.scss` must be manually compiled/ported to the corresponding `.css` file.

Site-wide configuration (menu, taxonomies, permalinks, author/social params, the homepage `whoami` terminal-block text, and `featuredProjects` slugs shown on the homepage) all lives in `hugo.toml`.
