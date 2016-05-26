NAME=cbcbSEQ
VERSION=0.9.1

export _R_CHECK_FORCE_SUGGESTS_=FALSE

all: clean document reference check build install test

install:
	echo "Performing R CMD INSTALL ${NAME}"
	ls && cd ../ && R CMD INSTALL ${NAME} && cd ${NAME}

reference:
	echo "Generating reference manual with R CMD Rd2pdf"
	rm -f inst/doc/reference.pdf
	R CMD Rd2pdf . -o inst/doc/reference.pdf

check:
	echo "Performing check with R CMD check ${NAME}"
	cd ../ && export _R_CHECK_FORCE_SUGGESTS_=FALSE && R CMD check ${NAME} --no-build-vignettes && cd ${NAME}

build:
	echo "Performing build with R CMD build ${NAME}"
	cd ../ && R CMD build ${NAME} && cd ${NAME}

test:
	echo "Running run_tests.R"
	./run_tests.R

roxygen:
	echo "Generating documentation with roxygen2::roxygenize()"
	Rscript -e "roxygen2::roxygenize()"

document:
	echo "Generating documentation with devtools::document()"
	Rscript -e "devtools::document()"

vignette:
	echo "Building vignettes with devtools::build_vignettes()"
	Rscript -e "devtools::build_vignettes()"

clean_vignette:
	rm -f inst/doc/*

vt:	clean_vignette vignette install

clean:
	rm -rf ${NAME}/
	rm -rf ${NAME}.Rcheck/
	rm -rf ${NAME}_${VERSION}.tar.gz
	find . -type f -name '*.Rdata' -exec rm -rf {} ';' 2>/dev/null

