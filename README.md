# Vulgaris — Alpha

A documentation-first portfolio built on [Hugo](https://gohugo.io), named for
*Phaseolus vulgaris* (the common bean). Terminal-tech dark theme, two content
types (`projects/`, `blog/`), tag-based "skills search," and a GitHub Actions
workflow that deploys straight to GitHub Pages.

## What's in here

```
content/
  projects/    → one Markdown file per repo/pipeline (didactic write-ups)
  blog/        → short log entries, notes, travelogues
layouts/       → all templates (home, list, single, taxonomy, 404)
static/        → css/style.css, js/main.js, img/favicon.svg
archetypes/    → front-matter templates for `hugo new`
.github/
  workflows/hugo.yml → builds + deploys to Pages on every push to main
```

