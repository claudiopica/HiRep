# Internal variables.
SPHINXOPTS   = -d _build/doctrees .

.PHONY: help
help:
	@echo "  clean       to remove generated documenation files."
	@echo "  html        to make standalone HTML files."
	@echo "  pdf         to make LaTeX files and run them through pdflatex."

.PHONY: clean
clean:
	rm -rf _build/* doxygen/build

.PHONY: html
html:
	sphinx-build -b html $(SPHINXOPTS) _build/html

.PHONY: pdf
pdf:
	sphinx-build -b latex $(SPHINXOPTS) _build/latex
	make -C _build/latex all-pdf
