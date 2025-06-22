root = .
include ${root}/defs.mk

pyprogs = $(shell file -F $$'\t' bin/* doc/bin/* | awk '/Python script/{print $$1}')

all:  lint

lint:
	${FLAKE8} --color=never ${pyprogs} lib/uniprotmap

tags: etags

etags:
	ctags -e lib/uniprotmap/*.py
test:
	cd tests && ${MAKE} test

fulltest:
	cd tests && ${MAKE} fulltest

clean:
	cd tests && ${MAKE} clean
	rm -rf  lib/uniprotmap/__pycache__ TAGS

savebak:
	savebak -git uniprotFeatMap
