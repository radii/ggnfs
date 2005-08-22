P1=ggnfs-doc

.PHONY: all clean squeaky

all : $(P1).pdf

$(P1).pdf : $(P1).tex
# Apparently, latex needs to be run at least once to create
# references and table of contents, which pdflatex seems to be having
# trouble with.
#	latex $(P1).tex
	pdflatex $(P1).tex
	pdflatex $(P1).tex
	mv $(P1).pdf ..

$(P1).ps : $(P1).tex
	convert reltst150.pdf reltst150.eps
	latex $(P1).tex
#	bibtex $(P1)
	latex $(P1).tex
	dvips -o $(P1).ps $(P1).dvi
	cp $(P1).ps .. ;
	rm reltst150.eps

clean :
	-rm -f *.aux *.log *.blg

squeaky : clean
	-rm -f $(P1).ps $(P1).dvi $(P1).bbl $(P1).pdf $(P1).toc reltst150.eps

