CORPORA_VERSION := $(shell perl -ne 'print $$1 if /^Version:\s+([0-9]+\.[0-9]+(-[0-9]+)?)$$/' corpora/DESCRIPTION)
CORPORA := corpora_$(CORPORA_VERSION).tar.gz

SIGIL_VERSION := $(shell perl -ne 'print $$1 if /^Version:\s+([0-9]+\.[0-9]+(-[0-9]+)?)$$/' SIGIL/DESCRIPTION)
SIGIL := SIGIL_$(SIGIL_VERSION).tar.gz

PKGS = $(CORPORA) $(SIGIL)

default:
	@echo "Available targets:"
	@echo "  make build[_pkg]    ...  build source package(s)"
	@echo "  make check[_pkg]    ...  test source package(s)"
	@echo "  make install[_pkg]  ...  install source package(s) locally"
	@echo "  make release[_pkg]  ...  check out and test source package(s) for CRAN release"
	@echo "  make binary[_pkg]   ...  compile binary package(s) for local architecture"
	@echo "  make clean[_pkg]    ...  clean up generated files"
	@echo "Packages: corpora (v$(CORPORA_VERSION)), SIGIL (v$(SIGIL_VERSION))"

.PHONY:	build check install release clean
.PHONY:	build_corpora check_corpora install_corpora release_corpora clean_corpora
.PHONY:	build_SIGIL check_SIGIL install_SIGIL release_SIGIL clean_SIGIL


build:	build_corpora build_SIGIL

build_corpora:
	R CMD build corpora

build_SIGIL:
	R CMD build SIGIL


check:	check_corpora check_SIGIL

check_corpora:	clean_corpora build_corpora $(CORPORA)
	R CMD check --multiarch --as-cran $(CORPORA)

check_SIGIL:	clean_SIGIL build_SIGIL $(SIGIL)
	R CMD check --multiarch --as-cran $(SIGIL)


install:	install_corpora install_SIGIL

install_corpora:	clean_corpora build_corpora $(CORPORA)
	R CMD INSTALL $(CORPORA)

install_SIGIL:	clean_SIGIL build_SIGIL $(SIGIL)
	R CMD INSTALL $(SIGIL)


release:	release_corpora release_SIGIL

release_corpora:	clean_corpora
	R CMD build --resave-data corpora
	R CMD check --as-cran $(CORPORA)
	@echo "$(CORPORA) is ready for release to CRAN"

release_SIGIL:	clean_SIGIL
	R CMD build --resave-data SIGIL
	R CMD check --as-cran $(SIGIL)
	@echo "$(SIGIL) is ready for release to CRAN"


binary:	binary_corpora binary_SIGIL

binary_corpora:	build_corpora $(CORPORA)
	R CMD INSTALL --build $(CORPORA)

binary_SIGIL:	build_SIGIL $(SIGIL)
	R CMD INSTALL --build $(SIGIL)


clean:	clean_corpora clean_SIGIL

clean_corpora:
	-rm -f corpora_*.tar.gz  # source package(s)
	-rm -f corpora_*.tgz     # Mac binaries
	-rm -f corpora_*.zip     # Windows binaries
	-rm -rf corpora.Rcheck   # working directory of R CMD check

clean_SIGIL:
	-rm -f SIGIL_*.tar.gz  # source package(s)
	-rm -f SIGIL_*.tgz     # Mac binaries
	-rm -f SIGIL_*.zip     # Windows binaries
	-rm -rf SIGIL.Rcheck   # working directory of R CMD check
