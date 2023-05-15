SHELL = /bin/bash -beEu -o pipefail

MAKEFLAGS += -rR
.PRECIOUS: %

algos = blat blast
swissprotGencodeFilts = raw paired
isoCanonFilts = raw paired

sppre = ../data/swissprot.9606
spsizes = ${sppre}.sizes
spmeta= ${sppre}.tab
spannot= ${sppre}.annots.tab
spfa = ${sppre}.fa.gz

gnpre = ../data/gencode.v43.pc_transcripts
gnsizes = ${gnpre}.sizes
gnmeta = ${gnpre}.meta.tsv
gnfa = ${gnpre}.fa.gz
gnpsl = ${gnpre}.psl


all: ${algos:%=%_mkAlgo}
	${MAKE} combineStatsCnts
%_mkAlgo:
	${MAKE} mkFilt algo=$*
	${MAKE} mkUniprotMapAnnots mkDecorators algo=$*

mkFilt: ${swissprotGencodeFilts:%=%_mkProtGeneFilt} ${isoCanonFilts:%=%_mkIsoCanonFilt}

%_mkProtGeneFilt: mkProtGeneAln
	${MAKE} mkProtGeneStats algo=${algo} filt=$*
%_mkIsoCanonFilt: mkIsoCanonAlgoAln
	${MAKE} mkIsoCanonStats algo=${algo} filt=$*


combineStatsCnts: swissprotGencode/swissprotGencode.stats-cnt.tsv \
		  isoCanon/isoCanon.stats-cnt.tsv

##
# swissprotGencode stats
# algo=, filt=
##
protGeneAlgoPre = swissprotGencode/swissprotGencode.${algo}
protGenePre = ${protGeneAlgoPre}.${filt}
protGeneAln = ${protGenePre}.psl

mkProtGeneStats: ${protGenePre}.qstats.tsv  ${protGenePre}.ostats.tsv ${protGenePre}.tstats.tsv \
	${protGenePre}.qstats-cnt.tsv ${protGenePre}.tstats-cnt.tsv

${protGenePre}.qstats.tsv: ${protGeneAln}
	pslStats -tsv -queryStats -queries=${spsizes} ${protGeneAln} $@

${protGenePre}.tstats.tsv: ${protGeneAln}
	pslSwap ${protGeneAln} stdout | pslStats -tsv -queryStats -queries=${gnsizes} /dev/stdin $@

${protGenePre}.ostats.tsv: ${protGeneAln}
	pslStats -tsv -overallStats -queries=${spsizes} ${protGeneAln} $@

${protGenePre}.%-cnt.tsv: ${protGenePre}.%.tsv
	./bin/countQStats swissprotGencode ${algo} ${filt} $* $< $@

# done in a second pass
swissprotGencode/swissprotGencode.stats-cnt.tsv: \
		$(foreach algo,${algos},$(foreach filt,${swissprotGencodeFilts},swissprotGencode/swissprotGencode.${algo}.${filt}.qstats-cnt.tsv)) \
		$(foreach algo,${algos},$(foreach filt,${swissprotGencodeFilts},swissprotGencode/swissprotGencode.${algo}.${filt}.tstats-cnt.tsv))
	tmlr cat $^ | tmlr sort -f stats -r filt > $@


##
# swissprotGencode alignment
#   algo=
##
mkProtGeneAln: ${protGeneAlgoPre}.raw.psl
protGeneAlnTmp = ${protGeneAlgoPre}.tmp
${protGeneAlgoPre}.raw.psl:
	@mkdir -p $(dir $@)
	rm -rf ${protGeneAlnTmp}
	../bin/proteinTranscriptAlign --algo=${algo} ${spfa} ${gnfa} $@.tmp ${protGeneAlnTmp}
	find ${protGeneAlnTmp}/aligns/ -name '*.psl' | xargs cat > ${protGeneAlgoPre}.realraw.psl
	mv -f $@.tmp $@

${protGeneAlgoPre}.paired.psl: ${protGeneAlgoPre}.raw.psl
	../bin/uniprotGencodeSelect ${gnmeta} ${spmeta} $< $@.tmp ${protGeneAlgoPre}.paired.problems.tsv
	pslCheck $@.tmp
	mv -f $@.tmp $@

##
# isoCanon stats
# algo=, filt=
## 
isoCanonAlgoPre = isoCanon/isoCanon.${algo}
isoCanonPre = ${isoCanonAlgoPre}.${filt}
isoCanonAln = ${isoCanonPre}.psl

mkIsoCanonStats: ${isoCanonPre}.qstats.tsv ${isoCanonPre}.ostats.tsv ${isoCanonPre}.qstats-cnt.tsv

${isoCanonPre}.qstats.tsv: ${isoCanonAln}
	pslStats -tsv -queryStats -queries=${spsizes} ${isoCanonAln} $@

${isoCanonPre}.ostats.tsv: ${isoCanonAln}
	pslStats -tsv -overallStats -queries=${spsizes} ${isoCanonAln} $@

${isoCanonPre}.%-cnt.tsv: ${isoCanonPre}.%.tsv
	./bin/countQStats isoCanon ${algo} ${filt} $* $< $@

# done in a second pass
isoCanon/isoCanon.stats-cnt.tsv: \
		$(foreach algo,${algos},$(foreach filt,${isoCanonFilts},isoCanon/isoCanon.${algo}.${filt}.qstats-cnt.tsv))
	tmlr cat $^ | tmlr sort -f stats -r filt > $@

##
# isoCanon alignment
#   algo=
##
mkIsoCanonAlgoAln: ${isoCanonAlgoPre}.raw.psl
isoCanonAlnTmp = ${isoCanonAlgoPre}.tmp

${isoCanonAlgoPre}.raw.psl:
	@mkdir -p $(dir $@)
	rm -rf ${isoCanonAlnTmp}
	../bin/uniprotIsoCanonicalAlign --algo=${algo} ${spfa} $@.tmp ${isoCanonAlnTmp}
	find ${isoCanonAlnTmp}/aligns/ -name '*.psl' | xargs cat > ${isoCanonAlgoPre}.realraw.psl
	mv -f $@.tmp $@

${isoCanonAlgoPre}.paired.psl: ${isoCanonAlgoPre}.raw.psl
	../bin/uniprotIsoCanonicalSelect ${spmeta} $< $@.tmp ${isoCanonAlgoPre}.paired.problems.tsv
	pslCheck $@.tmp
	mv -f $@.tmp $@

##
# map annotations
##
mapAnnotsAlgoPre = mapAnnots/mapAnnots.${algo}
mapAnnotsPsl = ${mapAnnotsAlgoPre}.psl
mapAnnotsRevTsv = ${mapAnnotsAlgoPre}.ref.tsv

mkUniprotMapAnnots: ${mapAnnotsPsl}

${mapAnnotsPsl}: ${spannot} ${isoCanonAlgoPre}.paired.psl ${protGeneAlgoPre}.paired.psl
	@mkdir -p $(dir $@)
	../bin/uniprotMapAnnots --keep ${spannot} ${isoCanonAlgoPre}.paired.psl ${protGeneAlgoPre}.paired.psl	${gnpsl} \
            $@.tmp ${mapAnnotsRevTsv}
	pslCheck $@.tmp
	mv -f $@.tmp $@

##
# make decorators
##
decoratorsAlgoPre = decorators/swissprot-gencode.v43.${algo}
decoratorsBed = ${decoratorsAlgoPre}.decorators.bed

mkDecorators: ${decoratorsBed}

${decoratorsBed}: ${spannot} ${mapAnnotsPsl}
	@mkdir -p $(dir $@)
	../bin/uniprotAnnotsToDecorators ${spannot} ${mapAnnotsPsl} ${mapAnnotsRevTsv} $@.tmp
	mv -f $@.tmp $@
