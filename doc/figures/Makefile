.PHONY: all

TARGETS = $(addsuffix .svg, $(basename $(shell ls *.mp)))

all: $(TARGETS)

%.svg: %.mp
	mpost -s 'outputformat="svg"' -s 'outputtemplate="%j-%c.%o"' $<
	mv $(addsuffix -ps.svg, $(basename $<)) $(addsuffix .svg, $(basename $<))
