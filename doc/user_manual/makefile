
TEXFLAGS  = 
TEXFILES := $(wildcard *.tex)

all: six.pdf

clean:
	rm -f six.aux six.fdb_latexmk six.fls six.idx six.ilg six.ind six.log six.lot six.out six.toc six.pdf six.synctex.gz

six.pdf: $(TEXFILES)
	pdflatex $(TEXFLAGS) six.tex
	pdflatex $(TEXFLAGS) six.tex
	pdflatex $(TEXFLAGS) six.tex
