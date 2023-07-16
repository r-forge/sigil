CORPORA_VERSION := $(shell perl -ne 'print $$1 if /^Version:\s+([0-9]+\.[0-9]+(-[0-9]+)?)$$/' corpora/DESCRIPTION)
CORPORA := corpora_$(CORPORA_VERSION).tar.gz

PKGS = $(CORPORA)

default:
	@echo "Available targets:"
	@echo "  make build[_pkg]    ...  build source package(s)"
	@echo "  make check[_pkg]    ...  test source package(s)"
	@echo "  make install[_pkg]  ...  install source package(s) locally"
	@echo "  make release[_pkg]  ...  check out and test source package(s) for CRAN release"
	@echo "  make binary[_pkg]   ...  compile binary package(s) for local architecture"
	@echo "  make clean[_pkg]    ...  clean up generated files"
	@echo "Packages: corpora (v$(CORPORA_VERSION))"

.PHONY:	build check install release clean
.PHONY:	build_corpora check_corpora install_corpora release_corpora clean_corpora

build:	build_corpora

build_corpora:
	R CMD build corpora


check:	check_corpora

check_corpora:	clean_corpora build_corpora $(CORPORA)
	R CMD check --multiarch --as-cran $(CORPORA)


install:	install_corpora

install_corpora:	clean_corpora build_corpora $(CORPORA)
	R CMD INSTALL $(CORPORA)


release:	release_corpora

release_corpora:	clean_corpora
	R CMD build --resave-data corpora
	R CMD check --as-cran $(CORPORA)
	@echo "$(CORPORA) is ready for release to CRAN"


binary:	binary_corpora

binary_corpora:	build_corpora $(CORPORA)
	R CMD INSTALL --build $(CORPORA)


clean:	clean_corpora

clean_corpora:
	-rm -f corpora_*.tar.gz  # source package(s)
	-rm -f corpora_*.tgz     # Mac binaries
	-rm -f corpora_*.zip     # Windows binaries
	-rm -rf corpora.Rcheck   # working directory of R CMD check
