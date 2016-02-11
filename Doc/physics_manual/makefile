

all: sixphys.pdf


sixphys.pdf: sixphys.tex refframe.tex
	pdflatex sixphys.tex
	bibtex sixphys.aux
	pdflatex sixphys.tex
	pdflatex sixphys.tex

sixphys-diff.pdf: sixphys.tex
	git show HEAD^:sixphys.tex > sixphys-old.tex
	latexdiff sixphys-old.tex sixphys.tex> sixphys-diff.tex
	pdflatex sixphys-diff.tex
	bibtex sixphys-diff.aux
	pdflatex sixphys-diff.tex
	pdflatex sixphys-diff.tex

sixphys.html: sixphys.tex refframe.tex
	htlatex sixphys.tex
