TARGET   = report
TEXFILES = $(wildcard *.tex)

.PRECIOUS: %.pdf
.PHONY: all pdf

all: $(TARGET).pdf

pdf: $(TARGET).pdf

$(TARGET).dvi: $(TARGET).tex $(TEXFILES)
	latex $<
	bibtex $(TARGET)
	latex $<
	latex $<

%.ps: %.dvi
	dvips -t letter -Ppdf -G0 $< -o $@

%.pdf: %.ps
	ps2pdf $<

clean:
	rm -rf *.bbl *.blg *.aux *.log *~ *.bak *.ps *.dvi *.log *.out *.tmp 

cleanall:
	rm -rf *.bbl *.blg *.aux *.log *~ *.bak *.ps *.dvi *.log *.pdf svnver.tex *.out *.tmp 

spell:
	ispell -f ispell.words -t *.tex
