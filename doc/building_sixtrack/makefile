
manName=building_sixtrack

all: $(manName).pdf

clean:
	rm -f $(manName).toc $(manName).out $(manName).aux $(manName).log $(manName).pdf

$(manName).pdf: $(manName).tex
	pdflatex $(manName).tex
	pdflatex $(manName).tex
