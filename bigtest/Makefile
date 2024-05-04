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

algos = blat
#blast
taxid = 9606

all: ${algos:%=%_mkAlgo}
	cd ../hub && ${MAKE}
%_mkAlgo:
	${MAKE} mkAlgo algo=$*
	${MAKE} mergeDecorators algo=$*

mkAlgo: ${uniprotDivisions:%=%_mkDivision}

%_mkDivision:
	${MAKE} mkDivision algo=${algo} div=$*

##
# processing on an division and alignment algorithm 
#   algo=
#   div=
##

mkDivision: mkProtTransAln mkProtTransMap mkUniprotMapAnnots mkDecorators

##
# protein/transcript alignment
#   algo=
#   div=
##
protTransAlnPre = protTransAln/${uniprotDataset}/protTransAln.${algo}
protTransRawAln = ${protTransAlnPre}.raw.psl
protTransAlnTmpDir = ${protTransAlnPre}.${tmpext}


mkProtTransAln: ${protTransRawAln} ${protTransPairedAln}
${protTransRawAln}: ${unitprotFa} ${gencodeFa}
	@mkdir -p $(dir $@)
	rm -rf ${protTransAlnTmpDir}
	${time} ${uniprotProteinTranscriptAlign} --algo=${algo} ${unitprotFa} ${gencodeFa} $@.${tmpext} ${protTransAlnTmpDir}
	mv -f $@.${tmpext} $@

##
# protein/transcript alignment
#   algo=
#   div=
##
protTransMapPre = protTransMap/${uniprotDataset}/protTransMap.${algo}
protTransPairedAln = ${protTransMapPre}.paired.psl
protTransPairProbTsv = ${protTransMapPre}.paired.problems.tsv

# uncomment to see debugging intermediates
# protTransInterDir = ${protTransMapPre}.inter
# protTransInterPre = ${protTransInterDir}/${algo}.inter.
# protTransInterOpt = --interPrefix=${protTransInterPre}

mkProtTransMap: ${protTransPairedAln}
${protTransPairedAln}: ${protTransRawAln}
	@mkdir -p $(dir $@) ${protTransInterDir}
	${time} ${uniprotProteinTranscriptMap} ${protTransInterOpt} ${gencodeMeta} ${gencodeGp} ${gencodePsl} ${uniprotMeta} $< $@.${tmpext} ${protTransPairProbTsv}
	pslCheck $@.${tmpext}
	mv -f $@.${tmpext} $@

##
# map annotations
##
mapAnnotsAlgoPre = mapAnnots/${uniprotDataset}/mapAnnots.${algo}
mapAnnotsPsl = ${mapAnnotsAlgoPre}.psl
mapAnnotsRefTsv = ${mapAnnotsAlgoPre}.ref.tsv
mapAnnotsProblemsTsv = ${mapAnnotsAlgoPre}.problems.tsv

# uncomment to see debugging intermediates
# mapAnnotsInterDir = ${mapAnnotsAlgoPre}.inter
# mapAnnotsInterPre = ${mapAnnotsInterDir}/${algo}.inter.
# mapAnnotsInterOpt = --interPrefix=${mapAnnotsInterPre}

mkUniprotMapAnnots: ${mapAnnotsPsl}

${mapAnnotsPsl}: ${protTransPairedAln} ${uniprotAnnot}
	@mkdir -p $(dir $@) ${mapAnnotsInterDir}
	${time} ${uniprotMapAnnots} ${uniprotAnnot} ${protTransPairedAln} ${gencodePsl} \
            $@.${tmpext} ${mapAnnotsRefTsv} ${mapAnnotsProblemsTsv} ${mapAnnotsInterOpt}
	pslCheck $@.${tmpext}
	mv -f $@.${tmpext} $@

##
# make decorators
##
decoratorsAlgoPre = decorators/${uniprotDataset}/uniprot-gencode.${gencodeVer}.${algo}
decoratorsBed = ${decoratorsAlgoPre}.decorators.bed
decoratorsTypesTsv = ${decoratorsAlgoPre}.decorators.types.tsv

mkDecorators: ${decoratorsBed}

${decoratorsBed}: ${upannot} ${mapAnnotsPsl}
	@mkdir -p $(dir $@)
	${time} ${uniprotAnnotsToDecorators} --dataset=${uniprotDataset} --nprocs=${nprocs} --batchSize=10000 --featTypesTsv=${decoratorsTypesTsv} \
	    ${uniprotMeta} ${uniprotAnnot} ${mapAnnotsPsl} ${mapAnnotsRefTsv} ${gencodePsl} $@.${tmpext}
	mv -f $@.${tmpext} $@

##
# mergeDecorators
# algo=
##
swissprotDecoratorsBed = decorators/SwissProt/uniprot-gencode.${gencodeVer}.${algo}.decorators.bed
tremblDecoratorsBed = decorators/TrEMBL/uniprot-gencode.${gencodeVer}.${algo}.decorators.bed

mergedDecoratorsBed = decorators/uniprot-gencode.${gencodeVer}.${algo}.decorators.bed


mergeDecorators: ${mergedDecoratorsBed}
${mergedDecoratorsBed}: ${swissprotDecoratorsBed} ${tremblDecoratorsBed}
	${uniprotDecoratorsMerge} --outBed=$@.${tmpext} $^
	mv -f $@.${tmpext} $@

# paranoid cleaning
clean:
	@echo "Error: to really clean up everything, use make realclean" >&2
	@exit 1

realclean:
	rm -rf protTransAln/ protTransMap/ mapAnnots/ decorators/
