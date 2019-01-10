
TOP := ..
DOC ?= .
BIBDIR := $(DOC)/bib
BIBSRC := $(wildcard $(BIBDIR)/*.bib)
PANDOC ?= pandoc
TOC    ?= $(DOC)/gh-md-toc

OUT    := $(TOP)/README.md
IN     := $(DOC)/README.md.in
BIBFIL := $(DOC)/README.bib
OUT_HTML := $(DOC)/README.html

DEPS = $(IN) $(DOC)/Makefile $(BIBFIL)

$(OUT): $(DEPS)
	@rm -f $@
	$(PANDOC) --bibliography $(BIBFIL) -t gfm-citations -s $< -o $@
	$(TOC) --insert $@
	@rm -f $@.orig.* $@.toc.*
	@chmod -w $@

all: $(OUT) $(OUT_HTML)

html: $(OUT_HTML)

$(OUT_HTML): $(DEPS)
	@rm -f $@
	$(PANDOC) --bibliography $(BIBFIL) -t html-citations -s $< -o $@
	@chmod -w $@

$(BIBFIL): $(BIBSRC)
	rm -f $@; cat $^ > $@ && chmod -w $@

clean:
	rm -f $(OUT) $(OUT_HTML) $(BIBFIL)

.PHONY: all html clean
