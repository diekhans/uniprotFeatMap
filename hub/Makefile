root = ..
include ${root}/defs.mk
include ${root}/xspecies.defs.mk

decoDir = ${root}/bigtest/decorators
bigGpAs =  ${HOME}/kent/src/hg/lib/big1GenePred.as

hub_url = https://hgwdev.gi.ucsc.edu/~markd/kznf/uniprotFeatMap/hub/hub.txt

algos = blat
algo = blat

all: ${algos:%=%_mkAlgo} mkGencodeBb mkDoc xspeciesHub
	${MAKE} check
%_mkAlgo:
	${MAKE} mkDecorators algo=$*

#
# algo=
decoratorsInBed = ${decoDir}/uniprot-gencode.${gencodeVer}.${algo}.decorators.bed
decoratorsOutBb = hg38/uniprot-gencode.${gencodeVer}.${algo}.bb


mkDecorators:  ${decoratorsOutBb}

${decoratorsOutBb}: ${decoratorsInBed}
	bedToBigBed -type=bed12+ -as=${uniprotDecoAs} -tab ${decoratorsInBed} ${chromSizes} $@.tmp
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
	mv -f $@.tmp $@

mkDoc: hg38/gencodeV45_blat_deco.html

hg38/gencodeV45_blat_deco.html: ${root}/doc/deco-desc.html
	@mkdir -p $(dir $@)
	cp -f $< $@

check:
	hubCheck ${hub_url}


##
# xspecies
##
xspeciesHub: ${asm_names:%=%_xspeciesHub}
%_xspeciesHub:
	${MAKE} xspeciesHubAsm  asmname=$* asmacc=$(call asmname_to_acc,$*)

##
# asmname
# asmacc
##
xhub_cat_bb = ${asmname}/${asmname}.cat.bb
xhub_deco_bb = ${asmname}/${asmname}.cat.blat.deco.bb
xhub_deco_html = ${asmname}/deco.html
xhub_trackDb = ${asmname}/trackDb.txt
xspeciesHubAsm: ${xhub_cat_bb} ${xhub_deco_bb} ${xhub_deco_html} ${xhub_trackDb}

${xhub_cat_bb}: ${xspecies_cat_bb}
	@mkdir -p $(dir $@)
	ln -f $< $@

${xhub_deco_bb}: ${xspecies_decorators_bed}
	@mkdir -p $(dir $@)
	bedToBigBed -type=bed12+ -as=${uniprotDecoAs} -tab $< ${xspecies_chromsizes} $@.tmp
	mv -f $@.tmp $@

${xhub_deco_html}: ${root}/doc/deco-desc.html
	@mkdir -p $(dir $@)
	cp -f $< $@

${xhub_trackDb}: common/xspecies.trackDb.txt
	@mkdir -p $(dir $@)
	sed -e 's/@genebb@/$(notdir ${xhub_cat_bb})/' -e 's/@decobb@/$(notdir ${xhub_deco_bb})/' $< > $@.tmp
	mv -f $@.tmp $@
