root = ..
include ${root}/defs.mk
SHELL = /bin/bash -beEu -o pipefail

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
##
## cross-species mapping
##
asm_names = \
    ponAbe-GCA_028885655.2 \
    ponPyg-GCA_028885625.2 \
    gorGor-GCA_029281585.2 \
    panPan-GCA_029289425.2 \
    panTro-GCA_028858775.2 \
    symSyn-GCA_028878055.2

asmname_to_acc = $(word 2,$(subst -, ,${1}))

xspeciesMap: ${asm_names:%=%_xspeciesMap}

%_xspeciesMap:
	${MAKE} xspeciesMapAsm  asmname=$* asmacc=$(call asmname_to_acc,$*) algo=blat div=swissprot

#
# asmname=
# asmacc=
#
src_asmname = hg38
src_genome = /hive/data/genomes/hg38/hg38.2bit
src_chromSizes = /hive/data/genomes/hg38/chrom.sizes

xspecies_url = https://public.gi.ucsc.edu/~pnhebbar/t2t_autosomes/assemblyHub
xspecies_rootdir=/hive/users/markd/gencode/projs/cls/cls-t2t-primate/hg38/data/assemblies

xspecies_catbb_url = ${xspecies_url}/${asmacc}/tm.bb
xspecies_genome = ${xspecies_rootdir}/${asmname}/${asmname}.2bit
xspecies_chromsizes_url = ${xspecies_url}/${asmacc}/chrom.sizes
src_xspecies_chains = ${xspecies_rootdir}/${asmname}/${src_asmname}-${asmname}.chain.gz

xspecies_dir = xspecies/${asmname}
xspecies_chromsizes = ${xspecies_dir}/${asmname}.sizes
xspecies_cat_pre = ${xspecies_dir}/${asmname}.cat
xspecies_cat_bb = ${xspecies_cat_pre}.bb
xspecies_cat_psl = ${xspecies_cat_pre}.psl
xspecies_cat_bed = ${xspecies_cat_pre}.bed
xspecies_cat_fa = ${xspecies_cat_pre}.fa

src_xspecies_rna_aln = ${xspecies_dir}/${src_asmname}-${asmname}.rna.psl
xspecies_mapannot_pre = ${xspecies_dir}/${src_asmname}-${asmname}.mapAnnots
xspecies_mapannot_psl = ${xspecies_mapannot_pre}.psl
xspecies_mapannot_ref_tsv = ${xspecies_mapannot_pre}.ref.tsv
xspecies_mapannot_problems_tsv = ${xspecies_mapannot_pre}.problems.tsv

xspecies_decorators_pre = ${xspecies_dir}/${src_asmname}-${asmname}.decorators
xspecies_decorators_bed = ${xspecies_decorators_pre}.bed


xspeciesMapAsm: ${xspecies_decorators_bed}

${xspecies_decorators_bed}: ${uniprotMeta} ${uniprotAnnot} ${xspecies_mapannot_psl} ${xspecies_mapannot_ref_tsv}
	${uniprotAnnotsToDecorators} --nprocs=${nprocs} --dataset=SwissProt ${uniprotMeta} ${uniprotAnnot} \
	    ${xspecies_mapannot_psl} ${xspecies_mapannot_ref_tsv} \
	    ${xspecies_cat_psl} $@.tmp
	mv -f $@.tmp $@

${xspecies_mapannot_psl} ${xspecies_mapannot_ref_tsv} ${xspecies_mapannot_problems_tsv} &: ${src_xspecies_rna_aln}
	${uniprotMapAnnots} --xspeciesTransPsl=${src_xspecies_rna_aln} \
	    ${uniprotAnnot} ${protTransPairedAln} ${src_xspecies_rna_aln}\
            ${xspecies_mapannot_psl}.tmp ${xspecies_mapannot_ref_tsv} ${xspecies_mapannot_problems_tsv}
	mv -f ${xspecies_mapannot_psl}.tmp ${xspecies_mapannot_psl}

${src_xspecies_rna_aln}: ${gencodePsl} ${gencodeFa} ${src_xspecies_chains} ${xspecies_cat_psl} ${xspecies_cat_fa} ${xspecies_cat_bigcat}
	@mkdir -p $(dir $@)
	${xspeciesTrans2TransMap} ${gencodePsl} ${gencodeFa} ${src_xspecies_chains} ${xspecies_cat_psl} ${xspecies_cat_fa} /dev/stdout | \
	    ${xspeciesGencode2CatFilter} /dev/stdin ${xspecies_cat_bb} $@.tmp
	mv -f $@.tmp $@

${xspecies_cat_fa}: ${xspecies_cat_bed} ${xspecies_genome}
	@mkdir -p $(dir $@)
	twoBitToFa -bed=${xspecies_cat_bed} ${xspecies_genome} $@.tmp
	mv -f $@.tmp $@

${xspecies_cat_bed}: ${xspecies_cat_psl}
	@mkdir -p $(dir $@)
	pslToBed ${xspecies_cat_psl} $@.tmp
	mv -f $@.tmp $@

${xspecies_chromsizes}:
	@mkdir -p $(dir $@)
	wget -nv -O $@.tmp ${xspecies_chromsizes_url}
	mv -f $@.tmp $@

${xspecies_cat_psl}: ${xspecies_cat_bb} ${xspecies_chromsizes}
	@mkdir -p $(dir $@)
	bigGenePredToGenePred ${xspecies_cat_bb} /dev/stdout | \
	  genePredToPsl ${xspecies_chromsizes} /dev/stdin $@.tmp
	mv -f $@.tmp $@

${xspecies_cat_bb}:
	@mkdir -p $(dir $@)
	bin/fixCatName ${xspecies_catbb_url} ${xspecies_genome} ${xspecies_cat_bb}.tmp
	mv -f $@.tmp $@

# paranoid cleaning
clean:
	@echo "Error: to really clean up everything, use make realclean" >&2
	@exit 1

realclean:
	rm -rf protTransAln/ protTransMap/ mapAnnots/ decorators/
