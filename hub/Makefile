SHELL = /bin/bash -beEu -o pipefail

MAKEFLAGS += -rR
.PRECIOUS: %

decoDir = ../bigtest/decorators
gencodePre = ../data//gencode.v43.
gencodeGp = ${gencodePre}.pc_transcripts.gp 
gencodeMeta = ${gencodePre}.pc_metadata.tsv
gencodeFa = ${gencodePre}.pc_transcripts.fa.gz
chromSizes = /hive/data/genomes/hg38/chrom.sizes
decoAs =  ../../kent/src/hg/lib/decoration.as

algos = blat blast


all: ${algos:%=%_mkAlgo}
%_mkAlgo:
	${MAKE} mkDecorators algo=$*

#
# algo=
decoratorsInBed = ../bigtest/decorators/swissprot-gencode.v43.${algo}.decorators.bed
decoratorsOutBb = hg38/swissprot-gencode.v43.${algo}.bb


mkDecorators:  ${decoratorsOutBb}

${decoratorsOutBb}: ${decoratorsInBed}
	bedToBigBed -type=bed12+ -as=${decoAs} -tab ${decoratorsInBed} ${chromSizes} $@.tmp
	mv -f $@.tmp $@


