# include with
#   root=../..
#   include ${root}/defs.mk

.PRECIOUS:
.SECONDARY:

SHELL = /bin/bash
export BASHOPTS = -beEu -o pipefail
MAKEFLAGS += -rR

PYTHON = python3
FLAKE8 = python3 -m flake8

export PYTHONPATH:=${rootDir}/lib:${PYTHONPATH}
export PYTHONWARNINGS=always

binDir = ${root}/bin
diff = diff -u

uniprotProteinTranscriptAlign = ${binDir}/uniprotProteinTranscriptAlign 
uniprotProteinTranscriptMap = ${binDir}/uniprotProteinTranscriptMap
uniprotAnnotsMap = ${binDir}/uniprotAnnotsMap
uniprotAnnotsToDecorators = ${binDir}/uniprotAnnotsToDecorators
uniprotDecoratorsMerge = ${binDir}/uniprotDecoratorsMerge
xspeciesTrans2TransMap = ${binDir}/xspeciesTrans2TransMap
xspeciesGencode2CatFilter = ${binDir}/xspeciesGencode2CatFilter
uniprotInfo = ${binDir}/uniprotInfo
interproProteinTranscriptAlign = ${binDir}/interproProteinTranscriptAlign
interproAnnotsMap = ${binDir}/interproAnnotsMap
interproAnnotsToDecorators = ${binDir}/interproAnnotsToDecorators

dataDir = ${root}/data

# GENCODE set
gencodeVer = v47
gencodePre = ${root}/data/gencode.${gencodeVer}
gencodeGp = ${gencodePre}.pc.gp 
gencodePsl = ${gencodePre}.pc.psl
gencodeMeta = ${gencodePre}.tsv
gencodeFa = ${gencodePre}.pc.fa

# UniProt: need div and taxid set for some variables
uniprotDivisions = swissprot trembl
uniprotDataset_swissprot = SwissProt
uniprotDataset_trembl = TrEMBL

uniprotDataset = ${uniprotDataset_${div}}
uniprotPrefix = ${dataDir}/${div}.${taxid}
uniprotMeta= ${uniprotPrefix}.tab
uniprotAnnot= ${uniprotPrefix}.annots.tab
unitprotFa = ${uniprotPrefix}.fa.gz

chromSizes = /hive/data/genomes/hg38/chrom.sizes
uniprotDecoAs =  ${root}/etc/uniprotDecoration.as
interproDecoAs =  ${root}/etc/interproDecoration.as
