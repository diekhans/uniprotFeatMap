SHELL = /bin/bash -beEu -o pipefail

MAKEFLAGS += -rR
.PRECIOUS: %

gnver = v44

decoDir = ../bigtest/decorators
gencodePre = ../data/gencode.${gnver}
gencodeGp = ${gencodePre}.pc_transcripts.gp 
gencodeMeta = ${gencodePre}.pc_transcripts.meta.tsv
gencodeFa = ${gencodePre}.pc_transcripts.fa.gz
chromSizes = /hive/data/genomes/hg38/chrom.sizes
decoAs =  ../etc/uniprotDecoration.as
bigGpAs =  ../../kent/src/hg/lib/bigGenePred.as


algos = blat blast


all: ${algos:%=%_mkAlgo} mkGencodeBb
%_mkAlgo:
	${MAKE} mkDecorators algo=$*

#
# algo=
decoratorsInBed = ../bigtest/decorators/swissprot-gencode.${gnver}.${algo}.decorators.bed
decoratorsOutBb = hg38/swissprot-gencode.${gnver}.${algo}.bb


mkDecorators:  ${decoratorsOutBb}

${decoratorsOutBb}: ${decoratorsInBed}
	bedToBigBed -type=bed12+ -as=${decoAs} -tab ${decoratorsInBed} ${chromSizes} $@.tmp
	mv -f $@.tmp $@


gencodePre = ../data/gencode.${gnver}.pc_transcripts

gencodeTrackPre = hg38/gencode.${gnver}

mkGencodeBb: ${gencodeTrackPre}.bb

${gencodeTrackPre}.bb: ${gencodeTrackPre}.bed
	bedToBigBed -type=bed12+8 -tab -as=${bigGpAs} $< ${chromSizes} $@.tmp
	mv -f $@.tmp $@

${gencodeTrackPre}.bed: ${gencodePre}.gp
	tmlr cut -f transcriptId,geneName,geneId ${gencodePre}.meta.tsv > hg38/gencode.names
	tmlr cut -f transcriptId,transcriptType ${gencodePre}.meta.tsv > hg38/gencode.types
	genePredToBigGenePred -geneNames=hg38/gencode.names -geneType=hg38/gencode.types $< /dev/stdout | sort -k1,1 -k2,2n > $@.tmp
	rm hg38/gencode.names hg38/gencode.types
	mv -f $@.tmp $@

