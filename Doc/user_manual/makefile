
manName=six

all: $(manName).pdf

clean:
	rm -f $(manName).blg $(manName).bbl $(manName).toc $(manName).out $(manName).aux $(manName).log $(manName).lot

$(manName).pdf: $(manName).tex
	pdflatex $(manName).tex
	pdflatex $(manName).tex
