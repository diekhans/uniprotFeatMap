# include with
#   root=../..
#   include ${root}/defs.mk

.PRECIOUS:

SHELL = /bin/bash
export BASHOPTS = -beEu -o pipefail

PYTHON = python3
FLAKE8 = python3 -m flake8

export PYTHONPATH:=${rootDir}/lib:${PYTHONPATH}
export PYTHONWARNINGS=always

binDir = ${root}/bin
diff = diff -u

uniprotProteinTranscriptAlign = ${binDir}/uniprotProteinTranscriptAlign 
uniprotProteinTranscriptMap = ${binDir}/uniprotProteinTranscriptMap
uniprotMapAnnots = ${binDir}/uniprotMapAnnots
uniprotAnnotsToDecorators = ${binDir}/uniprotAnnotsToDecorators

