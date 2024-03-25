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
taxid = 9606

up_divisions = swissprot trembl
up_dataset_swissprot = SwissProt
up_dataset_trembl = TrEMBL

gnver = v45
gnpre = ../data/gencode.${gnver}
gnmeta = ${gnpre}.meta.tsv
gnfa = ${gnpre}.pc_transcripts.fa.gz
gnpsl = ${gnpre}.psl
gngp = ${gnpre}.gp


all: ${algos:%=%_mkAlgo}
	cd ../hub && ${MAKE}
%_mkAlgo:
	${MAKE} mkAlgo algo=$*
	${MAKE} mergeDecorators algo=$*

mkAlgo: ${up_divisions:%=%_mkDivision}

%_mkDivision:
	${MAKE} mkDivision algo=${algo} div=$*

##
# processing on an division and alignment algorithm 
#   algo=
#   div=
##

uppre = ../data/${div}.${taxid}
upmeta= ${uppre}.tab
upannot= ${uppre}.annots.tab
upfa = ${uppre}.fa.gz
mkDivision: mkProtTransAln mkUniprotMapAnnots mkDecorators


##
# protein/transcript alignment
#   algo=
#   div=
##
dataset = ${up_dataset_${div}}
protTransPre = protTransAln/${dataset}/protTrans_Aln.${algo}
protTransRawAln = ${protTransPre}.raw.psl
protTransPairedAln = ${protTransPre}.paired.psl
protTransPairProbTsv = ${protTransPre}.paired.problems.tsv
protTransAlnTmpDir = ${protTransPre}.${tmpext}

# uncomment to see debugging intermediates
protTransInterDir = ${protTransAlgoPre}.inter
protTransInterPre = ${protTransInterDir}/${algo}.inter.
protTransInterOpt = --interPrefix=${protTransInterPre}


mkProtTransAln: ${protTransRawAln} ${protTransPairedAln}
${protTransRawAln}: ${upfa} ${gnfa}
	@mkdir -p $(dir $@) ${protTransInterDir}
	rm -rf ${protTransAlnTmpDir}
	${time} ${uniprotProteinTranscriptAlign} --algo=${algo} ${upfa} ${gnfa} $@.${tmpext} ${protTransAlnTmpDir}
	mv -f $@.${tmpext} $@

${protTransPairedAln}: ${protTransRawAln}
	${time} ${uniprotProteinTranscriptMap} ${gnmeta} ${gngp} ${gnpsl} ${upmeta} $< $@.${tmpext} ${protTransPairProbTsv}
	pslCheck $@.${tmpext}
	mv -f $@.${tmpext} $@

##
# map annotations
##
mapAnnotsAlgoPre = mapAnnots/${dataset}/mapAnnots.${algo}
mapAnnotsPsl = ${mapAnnotsAlgoPre}.psl
mapAnnotsRefTsv = ${mapAnnotsAlgoPre}.ref.tsv
mapAnnotsProblemsTsv = ${mapAnnotsAlgoPre}.problems.tsv

# uncomment to see debugging intermediates
mapAnnotsInterDir = ${mapAnnotsAlgoPre}.inter
mapAnnotsInterPre = ${mapAnnotsInterDir}/${algo}.inter.
mapAnnotsInterOpt = --interPrefix=${mapAnnotsInterPre}

mkUniprotMapAnnots: ${mapAnnotsPsl}

${mapAnnotsPsl}: ${upannot} ${protTransPairedAln}
	@mkdir -p $(dir $@) ${mapAnnotsInterDir}
	${time} ${uniprotMapAnnots} ${upannot} ${protTransPairedAln} ${gnpsl} \
            $@.${tmpext} ${mapAnnotsRefTsv} ${mapAnnotsProblemsTsv} ${mapAnnotsInterOpt}
	pslCheck $@.${tmpext}
	mv -f $@.${tmpext} $@

##
# make decorators
##
decoratorsAlgoPre = decorators/${dataset}/uniprot-gencode.${gnver}.${algo}
decoratorsBed = ${decoratorsAlgoPre}.decorators.bed
decoratorsTypesTsv = ${decoratorsAlgoPre}.decorators.types.tsv

mkDecorators: ${decoratorsBed}

${decoratorsBed}: ${upannot} ${mapAnnotsPsl}
	@mkdir -p $(dir $@)
	${time} ${uniprotAnnotsToDecorators} --dataset=${dataset} --nprocs=${nprocs} --batchSize=10000 --featTypesTsv=${decoratorsTypesTsv} \
	    ${upmeta} ${upannot} ${mapAnnotsPsl} ${mapAnnotsRefTsv} ${gnpsl} $@.${tmpext}
	mv -f $@.${tmpext} $@

##
# mergeDecorators
# algo=
##
swissprotDecoratorsBed = decorators/SwissProt/uniprot-gencode.${gnver}.${algo}.decorators.bed
tremblDecoratorsBed = decorators/TrEMBL/uniprot-gencode.${gnver}.${algo}.decorators.bed

mergedDecoratorsBed = decorators/uniprot-gencode.${gnver}.${algo}.decorators.bed


mergeDecorators: ${mergedDecoratorsBed}
${mergedDecoratorsBed}: ${swissprotDecoratorsBed} ${tremblDecoratorsBed}
	${uniprotDecoratorsMerge} --outBed=$@.${tmpext} $^
	mv -f $@.${tmpext} $@

# paranoid cleaning
clean:
	@echo "Error: to really clean up everything, use make realclean" >&2
	@exit 1

realclean:
	rm -rf prot_trans_aln mapAnnots decorators
