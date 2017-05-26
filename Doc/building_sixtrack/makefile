
manName=building_sixtrack

all: $(manName).pdf

clean:
	rm -f $(manName).toc $(manName).out $(manName).aux $(manName).log

manual: $(manName).tex
	pdflatex $(manName).tex
	pdflatex $(manName).tex
