all: FlukaIO.html FlukaIO.pdf

%.html: %.asciidoc
	asciidoc -b html5 -a theme=volnitsky -a pygments $<

%.pdf: %.asciidoc
	a2x -v -fpdf --fop $<

