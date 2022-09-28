# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
SPHINXOPTS   = -d _build/doctrees $(PAPEROPT_$(PAPER)) .

.PHONY: clean
clean:
	rm -rf _build/*

.PHONY: html
html:
	sphinx-build -b html $(SPHINXOPTS) _build/html

.PHONY: pdf
pdf:
	sphinx-build -b latex $(SPHINXOPTS) _build/latex
	make -C _build/latex all-pdf
