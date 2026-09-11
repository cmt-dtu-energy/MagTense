# MagTense TechManual

The manual is written in [reStructuredText](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html)
under [`source/`](source) and built with [Sphinx](https://www.sphinx-doc.org).

**The website is the committed HTML in this folder.**
<https://cmt-dtu-energy.github.io/MagTense> is served by GitHub Pages directly
from `docs/` on `master`, so editing the sources is not enough — the generated
HTML has to be rebuilt and committed. `.github/workflows/docs.yml` does that
automatically for anything that lands on `master`; the instructions below are
for building it yourself.

> The `pages-build-deployment` job in the Actions tab is GitHub's own Pages
> deploy step. It only publishes files that are already committed and does not
> run Sphinx.

## Prerequisites

Sphinx and the Read the Docs theme, in whichever environment you build from:

```bash
pip install -r docs/requirements.txt
```

They are not part of `magtense-env` by default. `conf.py` imports
`sphinx_rtd_theme` at module level, so an environment without it fails the
build with `No module named 'sphinx_rtd_theme'`.

## Building

From `docs/source`:

| Command | Effect |
| --- | --- |
| `make preview` / `make.bat preview` | Build into `source/_build/html` and print the path to open. Nothing outside `_build` is touched. |
| `make check` / `make.bat check` | Same, but with warnings as errors. Nothing is published. |
| `make publish` / `make.bat html` | Build, then replace the published HTML in `docs/`. |
| `make help` / `make.bat` | List all Sphinx targets. |

To just look at the pages, use `preview` — it is the only target that cannot
change anything under `docs/`.

On Windows use `make.bat`, on Linux and macOS use `make`. Both invoke Sphinx as
`python -m sphinx`, i.e. through the interpreter of the *active* environment,
rather than through whichever `sphinx-build` happens to be first on `PATH` —
picking up a different environment is the usual cause of a failed build. Set
the `SPHINXBUILD` environment variable to override.

Run `check` before `publish` if you want to see the warnings; `publish` builds
without `-W` so that a stray warning does not block a rebuild.

If `make.bat` misbehaves, Sphinx can always be called directly — this is
exactly what the targets above do, and it needs no `make` at all. From the
repository root:

```bat
python -m sphinx -b html docs\source docs\source\_build\html
```

Then open `docs\source\_build\html\index.html`.

Both scripts refuse to touch the published HTML unless the build succeeded.
They then delete the previous output before copying the new one in, so pages
whose source has been removed do not linger on the website. `README.md` and
`.nojekyll` are hand-maintained and are left alone.

## Publishing

`make publish` / `make.bat html` only updates your working tree. Commit the
result under `docs/` and get it onto `master` for it to reach the website:

```bash
git add docs && git commit -m "docs: rebuild the TechManual"
```

If you would rather let CI do it, just commit the changes under `docs/source/`
— the workflow rebuilds and commits the HTML when the change reaches `master`.

## Editing

- One page per topic under `source/`, referenced from the `toctree` in
  [`source/index.rst`](source/index.rst) or from a nested `toctree`.
- `conf.py` enables `sphinx.ext.autosectionlabel` **without**
  `autosectionlabel_prefix_document`, which means every section title has to be
  unique across the whole manual. Two pages with a section called `Python` will
  produce a duplicate-label warning and an ambiguous `:ref:`. This is why some
  headings are worded more specifically than they would otherwise need to be.
- Cross-reference a section by its title: ``:ref:`Micromagnetic output` ``.
- Any editor works; VS Code with the reStructuredText extension gives a live
  preview.
