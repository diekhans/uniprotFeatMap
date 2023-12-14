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

mkAlgo: mkProtGeneAln mkUniprotMapAnnots mkDecorators


##
# swissprotGencode alignment
#   algo=
##
protGenePre = swissprotGencode/swissprotGencode.${algo}
protGeneRawAln = ${protGenePre}.raw.psl
protGenePairedAln = ${protGenePre}.paired.psl
protGenePairProbTsv = ${protGenePre}.paired.problems.tsv
protGeneAlnTmpDir = ${protGenePre}.${tmpext}

mkProtGeneAln: ${protGeneRawAln} ${protGenePairedAln}
${protGeneRawAln}: ${spfa} ${gnfa}
	@mkdir -p $(dir $@)
	rm -rf ${protGeneAlnTmpDir}
	${time} ${uniprotProteinTranscriptAlign} --algo=${algo} ${spfa} ${gnfa} $@.${tmpext} ${protGeneAlnTmpDir}
	mv -f $@.${tmpext} $@

${protGenePairedAln}: ${protGeneRawAln}
	${time} ${uniprotProteinTranscriptMap} ${gnmeta} ${gngp} ${gnpsl} ${spmeta} $< $@.${tmpext} ${protGenePairProbTsv}
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
#mapAnnotsInterDir = ${mapAnnotsAlgoPre}.inter
#mapAnnotsInterPre = ${mapAnnotsInterDir}/${algo}.inter.
#mapAnnotsInterOpt = --interPrefix=${mapAnnotsInterPre}

mkUniprotMapAnnots: ${mapAnnotsPsl}

${mapAnnotsPsl}: ${spannot} ${protGenePairedAln}
	@mkdir -p $(dir $@) ${mapAnnotsInterDir}
	${time} ${uniprotMapAnnots} ${spannot} ${protGenePairedAln} ${gnpsl} \
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
