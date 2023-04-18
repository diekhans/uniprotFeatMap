.PRECIOUS:
SHELL = /bin/bash
export BASHOPTS = -beEu -o pipefail
PYTHON = python3
FLAKE8 = python3 -m flake8
pyprogs = $(shell file -F $$'\t' bin/* | awk '/Python script/{print $$1}')

all:  lint

lint:
	${FLAKE8} --color=never ${pyprogs}


