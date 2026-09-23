# Makefile for the MagTense TechManual.
#
#   make check     build with warnings as errors, without publishing
#   make publish   build and publish into ../ (the docs folder that GitHub
#                  Pages serves), then commit the result
#
# Any other target is passed straight through to Sphinx, e.g. "make latexpdf".
#
# Sphinx is invoked through the active interpreter rather than through whichever
# sphinx-build happens to come first on PATH: conf.py imports sphinx_rtd_theme,
# so picking up a different environment fails the build with
# "No module named 'sphinx_rtd_theme'". Override with
# "make publish SPHINXBUILD=sphinx-build" if you prefer the executable.

SPHINXOPTS    ?=
SPHINXBUILD   ?= python -m sphinx
SOURCEDIR     = .
BUILDDIR      = _build
PUBLISHDIR    = ..

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@echo
	@echo "MagTense targets:"
	@echo "  preview     build into _build/html for inspection, publishing nothing"
	@echo "  check       build with warnings as errors, publishing nothing"
	@echo "  publish     build and publish into the docs folder"

.PHONY: help preview check publish Makefile

# Build for inspection only. Nothing outside $(BUILDDIR) is touched, so this is
# always safe to run, and warnings do not stop it producing output.
preview:
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O)
	@echo
	@echo "Open this file in a browser:"
	@echo "    $(CURDIR)/$(BUILDDIR)/html/index.html"

check:
	@$(SPHINXBUILD) -b html -W -n "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build is clean. Nothing was published - run 'make publish' for that."

# The publish step runs only after a successful build; the old make.bat deleted
# the published assets before checking whether Sphinx had produced anything, so
# a failed build silently wiped the site.
publish:
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@test -f "$(BUILDDIR)/html/index.html" || \
		{ echo "No index.html was produced - $(PUBLISHDIR) left untouched."; exit 1; }
	@echo "Publishing into $(PUBLISHDIR) ..."
	@# Delete the previously published output first, so that pages whose source
	@# has been removed do not linger. README.md and .nojekyll are hand-maintained.
	@rm -rf "$(PUBLISHDIR)/_images" "$(PUBLISHDIR)/_sources" "$(PUBLISHDIR)/_static"
	@find "$(PUBLISHDIR)" -maxdepth 1 -type f \
		\( -name '*.html' -o -name '*.js' -o -name '*.inv' -o -name '.buildinfo' \) -delete
	@cp -a "$(BUILDDIR)/html/." "$(PUBLISHDIR)/"
	@rm -rf "$(BUILDDIR)"
	@echo "Done. Commit the changes under docs/ to update the website."

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
