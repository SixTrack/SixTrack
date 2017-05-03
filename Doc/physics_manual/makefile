
manName=sixphys
diffName=sixphys-diff

all: $(manName).pdf


clean:
	rm -f $(diffName).blg $(diffName).bbl $(diffName).toc $(diffName).out $(diffName).aux $(diffName).log $(diffName).tex
	rm -f sixphys-old.tex
	rm -f $(manName).blg $(manName).bbl $(manName).toc $(manName).out $(manName).aux $(manName).log 

cleanall:
	make clean
	rm -f $(diffName).pdf

manual: $(manName).tex refframe.tex
	pdflatex $(manName).tex
	bibtex $(manName).aux
	pdflatex $(manName).tex
	pdflatex $(manName).tex

diff: sixphys.tex
	git show HEAD^:Doc/physics_manual/sixphys.tex > sixphys-old.tex
	latexdiff sixphys-old.tex sixphys.tex> $(diffName).tex
	pdflatex $(diffName).tex
	bibtex $(diffName).aux
	pdflatex $(diffName).tex
	pdflatex $(diffName).tex

sixphys.html: sixphys.tex refframe.tex
	htlatex sixphys.tex
