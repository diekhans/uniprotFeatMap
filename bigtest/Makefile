SHELL = /bin/bash -beEu -o pipefail

MAKEFLAGS += -rR
.PRECIOUS: %

tmpext = $(shell hostname).$(shell echo $${PPID}).tmp


algos = blat blast
swissprotGencodeFilts = raw paired
isoCanonFilts = raw paired

sppre = ../data/swissprot.9606
spmeta= ${sppre}.tab
spannot= ${sppre}.annots.tab
spfa = ${sppre}.fa.gz

# rebuild this, as spfa changes
spsizes = ${sppre}.sizes

gnpre = ../data/gencode.v43.pc_transcripts
gnsizes = ${gnpre}.sizes
gnmeta = ${gnpre}.meta.tsv
gnfa = ${gnpre}.fa.gz
gnpsl = ${gnpre}.psl
nprocs = 64

all: ${algos:%=%_mkAlgo}
	${MAKE} combineStatsCnts
	cd ../hub && ${MAKE}
%_mkAlgo: ${spsizes}
	${MAKE} mkFilt algo=$*
	${MAKE} mkUniprotMapAnnots mkDecorators algo=$*

mkFilt: ${swissprotGencodeFilts:%=%_mkProtGeneFilt} ${isoCanonFilts:%=%_mkIsoCanonFilt}

%_mkProtGeneFilt: mkProtGeneAln
	${MAKE} mkProtGeneStats algo=${algo} filt=$*
%_mkIsoCanonFilt: mkIsoCanonAlgoAln
	${MAKE} mkIsoCanonStats algo=${algo} filt=$*

${spsizes}: ${spfa}
	faSize -detailed ${spfa} > $@.${tmpext}
	mv -f  $@.${tmpext} $@

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
protGeneAlnTmp = ${protGeneAlgoPre}.${tmpext}
${protGeneAlgoPre}.raw.psl:
	@mkdir -p $(dir $@)
	rm -rf ${protGeneAlnTmp}
	../bin/proteinTranscriptAlign --algo=${algo} ${spfa} ${gnfa} $@.${tmpext} ${protGeneAlnTmp}
	mv -f $@.${tmpext} $@

${protGeneAlgoPre}.paired.psl: ${protGeneAlgoPre}.raw.psl
	../bin/uniprotGencodeSelect ${gnmeta} ${spmeta} $< $@.${tmpext} ${protGeneAlgoPre}.paired.problems.tsv
	pslCheck $@.${tmpext}
	mv -f $@.${tmpext} $@

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
isoCanonAlnTmp = ${isoCanonAlgoPre}.${tmpext}

${isoCanonAlgoPre}.raw.psl:
	@mkdir -p $(dir $@)
	rm -rf ${isoCanonAlnTmp}
	../bin/uniprotIsoCanonicalAlign --algo=${algo} ${spfa} $@.${tmpext} ${isoCanonAlnTmp}
	find ${isoCanonAlnTmp}/aligns/ -name '*.psl' | xargs cat > ${isoCanonAlgoPre}.realraw.psl
	mv -f $@.${tmpext} $@

${isoCanonAlgoPre}.paired.psl: ${isoCanonAlgoPre}.raw.psl
	../bin/uniprotIsoCanonicalSelect ${spmeta} $< $@.${tmpext} ${isoCanonAlgoPre}.paired.problems.tsv
	pslCheck $@.${tmpext}
	mv -f $@.${tmpext} $@

##
# map annotations
##
mapAnnotsAlgoPre = mapAnnots/mapAnnots.${algo}
mapAnnotsPsl = ${mapAnnotsAlgoPre}.psl
mapAnnotsRevTsv = ${mapAnnotsAlgoPre}.ref.tsv

# uncomment to sae debugging intermediates
mapAnnotsInterDir = ${mapAnnotsAlgoPre}.inter
mapAnnotsInterPre = ${mapAnnotsInterDir}/${algo}.iter.
mapAnnotsInterOpt = --interPrefix=${mapAnnotsInterPre}

mkUniprotMapAnnots: ${mapAnnotsPsl}

${mapAnnotsPsl}: ${spannot} ${isoCanonAlgoPre}.paired.psl ${protGeneAlgoPre}.paired.psl
	@mkdir -p $(dir $@) ${mapAnnotsInterDir}
	../bin/uniprotMapAnnots ${spannot} ${isoCanonAlgoPre}.paired.psl ${protGeneAlgoPre}.paired.psl	${gnpsl} \
            $@.${tmpext} ${mapAnnotsRevTsv} ${mapAnnotsInterOpt}
	pslCheck $@.${tmpext}
	mv -f $@.${tmpext} $@

##
# make decorators
##
decoratorsAlgoPre = decorators/swissprot-gencode.v43.${algo}
decoratorsBed = ${decoratorsAlgoPre}.decorators.bed

mkDecorators: ${decoratorsBed}

${decoratorsBed}: ${spannot} ${mapAnnotsPsl}
	@mkdir -p $(dir $@)
	nice ../bin/uniprotAnnotsToDecorators --nprocs=${nprocs} --batchSize=10000 ${spannot} ${mapAnnotsPsl} ${mapAnnotsRevTsv} $@.${tmpext}
	mv -f $@.${tmpext} $@

# paranoid cleaning
clean:
	@echo "Error: to really clean up everything, use make realclean" >&2
	@exit 1

realclean:
	rm -rf swissprotGencode isoCanon mapAnnots decorators
