
TOP := ..
DOC ?= .
BIBDIR := $(DOC)/bib
BIBSRC := $(wildcard $(BIBDIR)/*.bib)
PANDOC ?= pandoc

OUT    := $(TOP)/README.md
IN     := $(DOC)/README.md.in
CSLFIL := $(DOC)/acm-siggraph.csl
BIBFIL := $(DOC)/README.bib

PANFLAGS = --bibliography $(BIBFIL) --csl $(CSLFIL) --filter pandoc-citeproc \
	   -r markdown -t markdown-citations -s $< -o $@

$(OUT): $(IN) $(BIBFIL) $(CSLFIL) $(DOC)/Makefile
	@rm -f $@
	$(PANDOC) $(PANFLAGS)
	@chmod -w $@

$(BIBFIL): $(BIBSRC)
	@rm -f $@ && cat $^ > $@ && chmod -w $@
