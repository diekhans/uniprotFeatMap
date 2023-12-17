root = ..
include ${root}/defs.mk

# used nprocs= to change number of processes
ifeq (${nprocs},)
  nprocs = 64
endif

MAKEFLAGS += -rR
.PRECIOUS: %

tmpext = $(shell hostname).$(shell echo $${PPID}).tmp

time = time -p

algos = blat blast

sppre = ../data/swissprot.9606
spmeta= ${sppre}.tab
spannot= ${sppre}.annots.tab
spfa = ${sppre}.fa.gz

gnver = v44
gnpre = ../data/gencode.${gnver}.pc_transcripts
gnmeta = ${gnpre}.meta.tsv
gnfa = ${gnpre}.fa.gz
gnpsl = ${gnpre}.psl
gngp = ${gnpre}.gp


all: ${algos:%=%_mkAlgo}
	cd ../hub && ${MAKE}
%_mkAlgo:
	${MAKE} mkAlgo algo=$*
	${MAKE} mkUniprotMapAnnots mkDecorators algo=$*

mkAlgo: mkProtTransAln mkUniprotMapAnnots mkDecorators


##
# swissprotGencode alignment
#   algo=
##
protTransPre = swissprotGencode/swissprotGencode.${algo}
protTransRawAln = ${protTransPre}.raw.psl
protTransPairedAln = ${protTransPre}.paired.psl
protTransPairProbTsv = ${protTransPre}.paired.problems.tsv
protTransAlnTmpDir = ${protTransPre}.${tmpext}

# uncomment to see debugging intermediates
protTransInterDir = ${protTransAlgoPre}.inter
protTransInterPre = ${protTransInterDir}/${algo}.inter.
protTransInterOpt = --interPrefix=${protTransInterPre}


mkProtTransAln: ${protTransRawAln} ${protTransPairedAln}
${protTransRawAln}: ${spfa} ${gnfa}
	@mkdir -p $(dir $@) ${protTransInterDir}
	rm -rf ${protTransAlnTmpDir}
	${time} ${uniprotProteinTranscriptAlign} --algo=${algo} ${spfa} ${gnfa} $@.${tmpext} ${protTransAlnTmpDir}
	mv -f $@.${tmpext} $@

${protTransPairedAln}: ${protTransRawAln}
	${time} ${uniprotProteinTranscriptMap} ${gnmeta} ${gngp} ${gnpsl} ${spmeta} $< $@.${tmpext} ${protTransPairProbTsv}
	pslCheck $@.${tmpext}
	mv -f $@.${tmpext} $@

##
# map annotations
##
mapAnnotsAlgoPre = mapAnnots/mapAnnots.${algo}
mapAnnotsPsl = ${mapAnnotsAlgoPre}.psl
mapAnnotsRefTsv = ${mapAnnotsAlgoPre}.ref.tsv
mapAnnotsProblemsTsv = ${mapAnnotsAlgoPre}.problems.tsv

# uncomment to see debugging intermediates
mapAnnotsInterDir = ${mapAnnotsAlgoPre}.inter
mapAnnotsInterPre = ${mapAnnotsInterDir}/${algo}.inter.
mapAnnotsInterOpt = --interPrefix=${mapAnnotsInterPre}

mkUniprotMapAnnots: ${mapAnnotsPsl}

${mapAnnotsPsl}: ${spannot} ${protTransPairedAln}
	@mkdir -p $(dir $@) ${mapAnnotsInterDir}
	${time} ${uniprotMapAnnots} ${spannot} ${protTransPairedAln} ${gnpsl} \
            $@.${tmpext} ${mapAnnotsRefTsv} ${mapAnnotsProblemsTsv} ${mapAnnotsInterOpt}
	pslCheck $@.${tmpext}
	mv -f $@.${tmpext} $@

##
# make decorators
##
decoratorsAlgoPre = decorators/swissprot-gencode.${gnver}.${algo}
decoratorsBed = ${decoratorsAlgoPre}.decorators.bed
decoratorsTypesTsv = ${decoratorsAlgoPre}.decorators.types.tsv

mkDecorators: ${decoratorsBed}

${decoratorsBed}: ${spannot} ${mapAnnotsPsl}
	@mkdir -p $(dir $@)
	${time} ${uniprotAnnotsToDecorators} --nprocs=${nprocs} --batchSize=10000 --featTypesTsv=${decoratorsTypesTsv} \
	    ${spmeta} ${spannot} ${mapAnnotsPsl} ${mapAnnotsRefTsv} ${gnpsl} $@.${tmpext}
	mv -f $@.${tmpext} $@

# paranoid cleaning
clean:
	@echo "Error: to really clean up everything, use make realclean" >&2
	@exit 1

realclean:
	rm -rf swissprotGencode mapAnnots decorators
