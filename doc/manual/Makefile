spatialmath.pdf: titlepage.pdf spatialmath.tex all.tex
	printf '\\def\\release{%s}\n' `cat ../../../RELEASE` > release.tex
	printf '\\def\\reldate{%s}\n' "`date '+%B %Y'`" >> release.tex
	pdflatex spatialmath
	sh funcidx_make.sh
	pdflatex spatialmath

titlepage.pdf: ../figs/tranimate.gif
	convert $<[10] $@

clean:
	-rm titlepage.pdf *.idx *.ind *.ilg *.aux *.log *.out *.toc

