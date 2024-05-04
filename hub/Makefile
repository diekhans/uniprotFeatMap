root = ..
include ${root}/defs.mk

decoDir = ../bigtest/decorators
bigGpAs =  ${HOME}/kent/src/hg/lib/bigGenePred.as

hub_url = https://hgwdev.gi.ucsc.edu/~markd/kznf/uniprotFeatMap/hub/hub.txt


all: ${algos:%=%_mkAlgo} mkGencodeBb mkDoc
	${MAKE} check
%_mkAlgo:
	${MAKE} mkDecorators algo=$*

#
# algo=
decoratorsInBed = ${decoDir}/uniprot-gencode.${gencodeVer}.${algo}.decorators.bed
decoratorsOutBb = hg38/uniprot-gencode.${gencodeVer}.${algo}.bb


mkDecorators:  ${decoratorsOutBb}

${decoratorsOutBb}: ${decoratorsInBed}
	bedToBigBed -type=bed12+ -as=${decoAs} -tab ${decoratorsInBed} ${chromSizes} $@.tmp
	mv -f $@.tmp $@


gencodeTrackPre = hg38/gencode.${gencodeVer}

mkGencodeBb: ${gencodeTrackPre}.bb

${gencodeTrackPre}.bb: ${gencodeTrackPre}.bed
	bedToBigBed -type=bed12+8 -tab -as=${bigGpAs} $< ${chromSizes} $@.tmp
	mv -f $@.tmp $@

${gencodeTrackPre}.bed: ${gencodeGp}
	tmlr cut -o -f transcriptId,geneName,geneId ${gencodeMeta} > hg38/gencode.names
	tmlr cut -o -f transcriptId,transcriptType ${gencodeMeta} > hg38/gencode.types
	genePredToBigGenePred -geneNames=hg38/gencode.names -geneType=hg38/gencode.types $< /dev/stdout | sort -k1,1 -k2,2n > $@.tmp
	#rm hg38/gencode.names hg38/gencode.types
	mv -f $@.tmp $@

mkDoc: hg38/gencodeV45_blat_deco.html

hg38/gencodeV45_blat_deco.html: ../doc/deco-desc.html
	cp -f $< $@

check: checked.done

#hg38/uniprot-gencode.${gencodeVer}.blast.bb
checked.done: ${gencodeTrackPre}.bb  hg38/uniprot-gencode.${gencodeVer}.blat.bb
	hubCheck ${hub_url}
	touch $@

