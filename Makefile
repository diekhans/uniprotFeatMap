root = .
include ${root}/defs.mk

pyprogs = $(shell file -F $$'\t' bin/* | awk '/Python script/{print $$1}')

all:  lint

lint:
	${FLAKE8} --color=never ${pyprogs} lib/uniprotmap

test:
	cd tests && ${MAKE} test

fulltest:
	cd tests && ${MAKE} fulltest

clean:
	cd tests && ${MAKE} clean
	rm -rf  lib/uniprotmap/__pycache__

savebak:
	savebak -git uniprotFeatMap
