DOTFLAGS=
DOT=dot

default: unfoc.pdf unfoc.png unfoc.svg

unfoc.pdf: unfoc.dot
	$(DOT) -Tpdf $? -o$@

unfoc.png: unfoc.dot
	$(DOT) -Tpng $? -o$@

unfoc.svg: unfoc.dot
	$(DOT) -Tsvg $? -o$@

clean:
	rm -f unfoc.pdf unfoc.png unfoc.svg
