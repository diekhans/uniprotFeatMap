SHELL = /bin/bash -beEu -o pipefail

MAKEFLAGS += -rR
.PRECIOUS: %

gnver = v45

decoDir = ../bigtest/decorators
decoAs =  ../etc/uniprotDecoration.as
bigGpAs =  ${HOME}/kent/src/hg/lib/bigGenePred.as

hub_url = https://hgwdev.gi.ucsc.edu/~markd/kznf/deco-hub/hub.txt

algos = blat
#blast


all: ${algos:%=%_mkAlgo} mkGencodeBb mkDoc
	${MAKE} check
%_mkAlgo:
	${MAKE} mkDecorators algo=$*

#
# algo=
decoratorsInBed = ../bigtest/decorators/uniprot-gencode.${gencodeVer}.${algo}.decorators.bed
decoratorsOutBb = hg38/uniprot-gencode.${gencodeVer}.${algo}.bb


mkDecorators:  ${decoratorsOutBb}

${decoratorsOutBb}: ${decoratorsInBed}
	bedToBigBed -type=bed12+ -as=${decoAs} -tab ${decoratorsInBed} ${chromSizes} $@.tmp
	mv -f $@.tmp $@


gencodePre = ../data/gencode.${gencodeVer}.pc_transcripts

gencodeTrackPre = hg38/gencode.${gencodeVer}

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

mkDoc: hg38/gencodeV45_blat_deco.html

hg38/gencodeV45_blat_deco.html: ../doc/deco-desc.html
	cp -f $< $@

check: checked.done

#hg38/uniprot-gencode.${gencodeVer}.blast.bb
checked.done: ${gencodeTrackPre}.bb  hg38/uniprot-gencode.${gencodeVer}.blat.bb
	hubCheck ${hub_url}
	touch $@
