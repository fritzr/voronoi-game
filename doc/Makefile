
TOP := ..
DOC ?= .
BIBDIR := $(DOC)/bib
BIBSRC := $(wildcard $(BIBDIR)/*.bib)
PANDOC ?= pandoc

OUT    := $(TOP)/README.md
IN     := $(DOC)/README.md.in
BIBFIL := $(DOC)/README.bib

PANFLAGS = --bibliography $(BIBFIL) --atx-headers -t gfm-citations -s $< -o $@

$(OUT): $(IN) $(BIBFIL) $(CSLFIL) $(DOC)/Makefile
	@rm -f $@
	$(PANDOC) $(PANFLAGS)
	@chmod -w $@

$(BIBFIL): $(BIBSRC)
	@rm -f $@ && cat $^ > $@ && chmod -w $@
