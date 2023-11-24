SHELL = /bin/bash -beEu -o pipefail

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
nprocs = 64

all: ${algos:%=%_mkAlgo}
	cd ../hub && ${MAKE}
%_mkAlgo:
	${MAKE} mkAlgo algo=$*
	${MAKE} mkUniprotMapAnnots mkDecorators algo=$*

mkAlgo: mkProtGeneAln mkUniprotMapAnnots mkDecorators


##
# swissprotGencode alignment
#   algo=
##
protGeneAlgoPre = swissprotGencode/swissprotGencode.${algo}
protGenePre = ${protGeneAlgoPre}.${filt}
protGeneAln = ${protGenePre}.psl

mkProtGeneAln: ${protGeneAlgoPre}.raw.psl
protGeneAlnTmp = ${protGeneAlgoPre}.${tmpext}
${protGeneAlgoPre}.raw.psl:
	@mkdir -p $(dir $@)
	rm -rf ${protGeneAlnTmp}
	${time} ../bin/uniprotProteinTranscriptAlign --algo=${algo} ${spfa} ${gnfa} $@.${tmpext} ${protGeneAlnTmp}
	mv -f $@.${tmpext} $@

${protGeneAlgoPre}.paired.psl: ${protGeneAlgoPre}.raw.psl
	${time} ../bin/uniprotGencodeSelect ${gnmeta} ${spmeta} $< $@.${tmpext} ${protGeneAlgoPre}.paired.problems.tsv
	pslCheck $@.${tmpext}
	mv -f $@.${tmpext} $@

##
# map annotations
##
mapAnnotsAlgoPre = mapAnnots/mapAnnots.${algo}
mapAnnotsPsl = ${mapAnnotsAlgoPre}.psl
mapAnnotsRefTsv = ${mapAnnotsAlgoPre}.ref.tsv

# uncomment to see debugging intermediates
#mapAnnotsInterDir = ${mapAnnotsAlgoPre}.inter
#mapAnnotsInterPre = ${mapAnnotsInterDir}/${algo}.inter.
#mapAnnotsInterOpt = --interPrefix=${mapAnnotsInterPre}

mkUniprotMapAnnots: ${mapAnnotsPsl}

${mapAnnotsPsl}: ${spannot} ${protGeneAlgoPre}.paired.psl
	@mkdir -p $(dir $@) ${mapAnnotsInterDir}
	${time} ../bin/uniprotMapAnnots ${spannot} ${protGeneAlgoPre}.paired.psl ${gnpsl} \
            $@.${tmpext} ${mapAnnotsRefTsv} ${mapAnnotsInterOpt}
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
	${time} ../bin/uniprotAnnotsToDecorators --nprocs=${nprocs} --batchSize=10000 --featTypesTsv=${decoratorsTypesTsv} \
	    ${spmeta} ${spannot} ${mapAnnotsPsl} ${mapAnnotsRefTsv} ${gnpsl} $@.${tmpext}
	mv -f $@.${tmpext} $@

# paranoid cleaning
clean:
	@echo "Error: to really clean up everything, use make realclean" >&2
	@exit 1

realclean:
	rm -rf swissprotGencode mapAnnots decorators
