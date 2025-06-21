root = ..
include ${root}/defs.mk

# see bin/getTestData for updating test data


# used nproc= to change number of processes
ifeq (${nproc},)
  nproc = 8
endif

all:
	@echo "Note: targets are:"
	@echo "      make test"
	@echo "      make fulltest"
	@echo "      make aligntest"
	@echo "      make hub"
	@echo "      make clean"
	@exit 1


test: uniprotTests interproTests xspeciesAnalyzeTests hub
	@echo "==========================================" >&2
	@echo "Note to test alignments use: make fulltest" >&2
	@echo "==========================================" >&2

fulltest: aligntest test

hg38ChromSizes = input/hg38.sizes

testHubFiles = empty.html genomes.txt hub.txt hg38.trackDb.txt \
	ponAbe-GCA_028885655.2.trackDb.txt \
	symSyn-GCA_028878055.2.trackDb.txt
testHubFilesOutput = ${testHubFiles:%=output/%}

testHubUrl = https://hgwdev.gi.ucsc.edu/~markd/kznf/uniprotFeatMap/tests/output/hub.txt

swissprotMeta = input/swissprot.9606.tab
swissprotAnnots = input/swissprot.9606.annots.tab

gencodeGp = input/gencode.v43.pc.gp
gencodePcFa = input/gencode.v43.pc.fa
gencodePcPsl = input/gencode.v43.pc.psl
gencodeProtFa = input/gencode.v43.proteins.fa
gencodeMeta = input/gencode.v43.metadata.tsv
gencodeBb = output/gencode.v43.bb
gencodeBrowserDir = output/gencode-dir
gencodeBrowserDirDone = ${gencodeBrowserDir}/done
gencodeInterproTsv = input/gencode.v43.interproscan.tsv

ponAbeCat1BigCat = input/ponAbe-GCA_028885655.2.cat.bigCat
ponAbeCat1Gp = input/ponAbe-GCA_028885655.2.cat.gp
ponAbeCat1Psl = input/ponAbe-GCA_028885655.2.cat.psl
ponAbeCat1RnaFa = input/ponAbe-GCA_028885655.2.cat.rna.fa
ponAbeCat1ProtFa = input/ponAbe-GCA_028885655.2.cat.prot.fa
ponAbeCat1Meta = input/ponAbe-GCA_028885655.2.cat.meta.tsv
ponAbeCat1Bb = output/ponAbe-GCA_028885655.2.cat.bb
ponAbeCat1InterproTsv = input/ponAbe-GCA_028885655.2.interpro.tsv
ponAbeTwobit = https://hgdownload.soe.ucsc.edu/hubs/GCA/028/885/655/GCA_028885655.2/GCA_028885655.2.2bit
ponAbeChains = input/ponAbe-GCA_028885655.2.chains
ponAbeChromSizes = input/ponAbe-GCA_028885655.2.sizes

symSynCat1BigCat = input/symSyn-GCA_028878055.2.cat.bigCat
symSynCat1Gp = input/symSyn-GCA_028878055.2.cat.gp
symSynCat1Psl = input/symSyn-GCA_028878055.2.cat.psl
symSynCat1RnaFa = input/symSyn-GCA_028878055.2.cat.rna.fa
symSynCat1ProtFa = input/symSyn-GCA_028878055.2.cat.prot.fa
symSynCat1Meta = input/symSyn-GCA_028878055.2.cat.meta.tsv
symSynCat1Bb = output/symSyn-GCA_028878055.2.cat.bb
symSynTwobit = https://hgdownload.soe.ucsc.edu/hubs/GCA/028/878/055/GCA_028878055.2/GCA_028878055.2.2bit
symSynChains = input/symSyn-GCA_028878055.2.chains
symSynChromSizes = input/symSyn-GCA_028878055.2.sizes

uniprotInterproAnnotAssocTsv = input/uniprotInterproAnnotAssoc.tsv


bigCatAs = input/bigCat.as

# alignment tests must be done serially due to parasol losing small jobs on multiple
# batches with short jobs
aligntest:
	${MAKE} protTransAlignTests
	${MAKE} interproAlignTests

####
uniprotTests: uniprotProteinTranscriptMapTests uniprotAnnotsMapTests uniprotAnnotsToDecoratorsTests \
	uniprotDecoratorsMergeTests xspeciesPonAbeTests xspeciesSymSynTests uniprotInfoTests 


# serial due to parasol bug
protTransAlignTests:
	${MAKE} testProtTransSPBlastAlign
	${MAKE} testProtTransSPBlatAlign
	${MAKE} testProtTransTRBlastAlign
	${MAKE} testProtTransTRBlatAlign

# $(call runUpAlignments,blat,swissprot)
define runUpAlignments
	rm -rf output/$@.psl output/$@.work
	${uniprotProteinTranscriptAlign} --algo=${1} input/${2}.9606.tab input/${2}.9606.fa ${gencodePcFa} output/$@.psl output/$@.work
	diff expected/$@.psl output/$@.psl
	pslCheck -verbose=0 output/$@.psl
endef

testProtTransSPBlastAlign: mkout
	$(call runUpAlignments,blast,swissprot)

testProtTransSPBlatAlign: mkout
	$(call runUpAlignments,blat,swissprot)

testProtTransTRBlastAlign: mkout
	$(call runUpAlignments,blast,trembl)

testProtTransTRBlatAlign: mkout
	$(call runUpAlignments,blat,trembl)

####
uniprotProteinTranscriptMapTests: testUniprotProteinTranscriptMapSP testUniprotProteinTranscriptMapTR

# turn on for debugging
uniprotProteinTranscriptMapInterPrefix = --interPrefix=output/$@.inter.
# uniprotProteinTranscriptMapInterPrefix =

# $call runUniprotProteinTranscriptMap,swissprot,SP)
define runUniprotProteinTranscriptMap
	${uniprotProteinTranscriptMap} ${uniprotProteinTranscriptMapInterPrefix} \
		GENCODE ${gencodeMeta} ${gencodePcPsl} ${gencodeGp} \
		input/${1}.9606.tab expected/testProtTrans${2}BlatAlign.psl \
		output/$@.psl output/$@.problems.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.problems.tsv output/$@.problems.tsv
	pslCheck -verbose=0 output/$@.psl
endef

testUniprotProteinTranscriptMapSP: mkout
	$(call runUniprotProteinTranscriptMap,swissprot,SP)

testUniprotProteinTranscriptMapTR: mkout
	$(call runUniprotProteinTranscriptMap,trembl,TR)

###
uniprotAnnotsMapTests: testUniprotAnnotsMapSP testUniprotAnnotsMapTR
# intermediates saved for debugging

# $(call runMapAnnots,swissprot,SP)
define runMapAnnots
	${uniprotAnnotsMap} --interPrefix=output/$@.inter. \
	    --annot2TransPsl=output/$@.annotsTrans.psl \
	    ${gencodePcPsl} \
	    input/${1}.9606.annots.tab expected/testUniprotProteinTranscriptMap${2}.psl \
            output/$@.psl output/$@.ref.tsv output/$@.problems.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.ref.tsv output/$@.ref.tsv
	diff expected/$@.problems.tsv output/$@.problems.tsv
	diff expected/$@.inter.annotTrans.psl output/$@.inter.annotTrans.psl
	diff expected/$@.inter.annotTrans.psl output/$@.annotsTrans.psl
	diff expected/$@.annotsTrans.psl output/$@.annotsTrans.psl
	pslCheck -verbose=0 output/$@.psl
	pslCheck -verbose=0 output/$@.inter.annotTrans.psl
endef

testUniprotAnnotsMapSP: mkout
	$(call runMapAnnots,swissprot,SP)
testUniprotAnnotsMapTR: mkout
	$(call runMapAnnots,trembl,TR)

###
uniprotAnnotsToDecoratorsTests: testUniprotAnnotsToDecoratorsSP testUniprotAnnotsToDecoratorsTR \
	testUniprotAnnotsToDecoColorHelp testUniprotAnnotsToDecoratorsError

# $(call runUniprotAnnotsToDecorators,swissprot,SP)
define runUniprotAnnotsToDecorators
	${uniprotAnnotsToDecorators} --nproc=${nproc} ${gencodePcPsl} \
	    input/$(1).9606.tab input/$(1).9606.annots.tab \
	    expected/testUniprotAnnotsMap$(2).psl expected/testUniprotAnnotsMap$(2).ref.tsv \
            output/$@.bed --featTypesTsv=output/$@.types.tsv
	diff expected/$@.bed output/$@.bed
	diff expected/$@.types.tsv output/$@.types.tsv
endef

testUniprotAnnotsToDecoratorsSP: mkout ${gencodeBb}
	$(call runUniprotAnnotsToDecorators,swissprot,SP)

testUniprotAnnotsToDecoratorsTR: mkout ${gencodeBb}
	$(call runUniprotAnnotsToDecorators,trembl,TR)

testUniprotAnnotsToDecoColorHelp: mkout
	${uniprotAnnotsToDecorators} --help-colors >output/$@.out
	diff expected/$@.out output/$@.out

# error handling for pools, use to not terminate on fail to open files
ifeq (${nproc},1)
testUniprotAnnotsToDecoratorsError:
else
testUniprotAnnotsToDecoratorsError: mkout ${gencodeBb}
	! ${uniprotAnnotsToDecorators} --nproc=${nproc} ${gencodePcPsl} no_swissprotMeta ${swissprotAnnots} \
	    expected/testUniprotAnnotsMapSP.psl expected/testUniprotAnnotsMapSP.ref.tsv \
            output/$@.bed >& output/$@.err || { echo "Error: command should fail"; exit 1; }
	grep -q 'Exception: Pool initialization failed' output/$@.err 
endif

###
uniprotDecoratorsMergeTests: testUniprotDecoratorsMerge

testUniprotDecoratorsMerge: mkout
	${uniprotDecoratorsMerge} --outBed=output/$@.bed expected/testUniprotAnnotsToDecoratorsSP.bed expected/testUniprotAnnotsToDecoratorsTR.bed
	bedToBigBed -type=bed12+ -as=${uniprotDecoAs} -tab output/$@.bed ${hg38ChromSizes} output/$@.bb
	diff expected/$@.bed output/$@.bed

###
# cross-species mapping with orangutan
###
xspeciesPonAbeTests: testXsPonAbeTransMap testXsPonAbeRnaMap testXsPonAbeGencodeCatFilter testXsPonAbeUniprotAnnotsMap \
	testXsPonAbeAnnotsToDeco

testXsPonAbeTransMap: mkout
	${xspeciesTrans2GenomeMap} ${gencodePcPsl} ${gencodePcFa} \
		${ponAbeChains} ${ponAbeTwobit} output/$@.psl
	diff expected/$@.psl output/$@.psl

testXsPonAbeRnaMap: mkout
	${xspeciesTrans2TransMap} expected/testXsPonAbeTransMap.psl ${gencodePcFa} \
		${ponAbeCat1Psl} ${ponAbeCat1RnaFa} output/$@.psl
	diff expected/$@.psl output/$@.psl

testXsPonAbeGencodeCatFilter: mkout ${ponAbeCat1Bb}
	${xspeciesGencode2CatFilter} expected/testXsPonAbeRnaMap.psl ${ponAbeCat1Bb} output/$@.psl
	diff expected/$@.psl output/$@.psl

testXsPonAbeUniprotAnnotsMap: mkout
	${uniprotAnnotsMap} --interPrefix=output/$@.inter. \
	    --annot2TransPsl=output/$@.annotsTrans.psl \
	    --xspeciesTrans2TransPsl=expected/testXsPonAbeGencodeCatFilter.psl \
	    ${ponAbeCat1Psl} ${swissprotAnnots} expected/testUniprotProteinTranscriptMapSP.psl \
            output/$@.psl output/$@.ref.tsv output/$@.problems.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.ref.tsv output/$@.ref.tsv
	diff expected/$@.problems.tsv output/$@.problems.tsv
	diff expected/$@.inter.annotTrans.psl output/$@.inter.annotTrans.psl
	diff expected/$@.inter.annotTrans.psl output/$@.annotsTrans.psl
	diff expected/$@.annotsTrans.psl output/$@.annotsTrans.psl
	pslCheck -verbose=0 output/$@.psl
	pslCheck -verbose=0 output/$@.inter.annotTrans.psl

testXsPonAbeAnnotsToDeco: mkout
	${uniprotAnnotsToDecorators} --nproc=${nproc} ${ponAbeCat1Psl} ${swissprotMeta} ${swissprotAnnots} \
	    expected/testXsPonAbeUniprotAnnotsMap.psl expected/testXsPonAbeUniprotAnnotsMap.ref.tsv \
            output/$@.bed --featTypesTsv=output/$@.types.tsv
	diff expected/$@.bed output/$@.bed
	diff expected/$@.types.tsv output/$@.types.tsv
	bedToBigBed -type=bed12+ -as=${uniprotDecoAs} -sizesIs2Bit -tab expected/$@.bed ${ponAbeTwobit} output/$@.bb

##
# Cross-species mapping with gibbon
# Gibbon has a domain deletion in ZNF649
#
xspeciesSymSynTests: testXsSymSynTransMap testXsSymSynRnaMap testXsSymSynGencodeCatFilter testXsSymSynUniprotAnnotsMap \
	testXsSymSynAnnotsToDeco

testXsSymSynTransMap: mkout
	${xspeciesTrans2GenomeMap} ${gencodePcPsl} ${gencodePcFa} \
		${symSynChains} ${symSynTwobit} output/$@.psl
	diff expected/$@.psl output/$@.psl

testXsSymSynRnaMap: mkout
	${xspeciesTrans2TransMap} expected/testXsSymSynTransMap.psl ${gencodePcFa} \
		${symSynCat1Psl} ${symSynCat1RnaFa} output/$@.psl
	diff expected/$@.psl output/$@.psl

testXsSymSynGencodeCatFilter: mkout ${symSynCat1Bb}
	${xspeciesGencode2CatFilter} expected/testXsSymSynRnaMap.psl ${symSynCat1Bb} output/$@.psl
	diff expected/$@.psl output/$@.psl

testXsSymSynUniprotAnnotsMap: mkout
	${uniprotAnnotsMap} --interPrefix=output/$@.inter. \
	    --annot2TransPsl=output/$@.annotsTrans.psl \
	    --xspeciesTrans2TransPsl=expected/testXsSymSynGencodeCatFilter.psl \
	    ${symSynCat1Psl} ${swissprotAnnots} expected/testUniprotProteinTranscriptMapSP.psl \
            output/$@.psl output/$@.ref.tsv output/$@.problems.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.ref.tsv output/$@.ref.tsv
	diff expected/$@.problems.tsv output/$@.problems.tsv
	diff expected/$@.inter.annotTrans.psl output/$@.inter.annotTrans.psl
	diff expected/$@.inter.annotTrans.psl output/$@.annotsTrans.psl
	diff expected/$@.annotsTrans.psl output/$@.annotsTrans.psl
	pslCheck -verbose=0 output/$@.psl
	pslCheck -verbose=0 output/$@.inter.annotTrans.psl

testXsSymSynAnnotsToDeco: mkout
	${uniprotAnnotsToDecorators} --nproc=${nproc} ${symSynCat1Psl} ${swissprotMeta} ${swissprotAnnots} \
	    expected/testXsSymSynUniprotAnnotsMap.psl expected/testXsSymSynUniprotAnnotsMap.ref.tsv \
            output/$@.bed --featTypesTsv=output/$@.types.tsv
	diff expected/$@.bed output/$@.bed
	diff expected/$@.types.tsv output/$@.types.tsv
	bedToBigBed -type=bed12+ -as=${uniprotDecoAs} -sizesIs2Bit -tab expected/$@.bed ${symSynTwobit} output/$@.bb


###
uniprotInfoTests: testUniprotInfoEnsTrans testUniprotInfoGeneName

testUniprotInfoEnsTrans: mkout
	${uniprotInfo} ensemblTransId ${swissprotMeta} output/$@.tsv
	diff expected/$@.tsv output/$@.tsv

testUniprotInfoGeneName:
	${uniprotInfo} geneName ${swissprotMeta} output/$@.tsv
	diff expected/$@.tsv output/$@.tsv

####
# serial due to parasol bug
interproAlignTests:
	${MAKE} testInterproBlatAlign
	${MAKE} testInterproPonAbeCat1BlatAlign

testInterproBlatAlign: mkout
	rm -rf output/$@.psl output/$@.work
	${interproProteinTranscriptAlign} GENCODE ${gencodeMeta} ${gencodePcFa} ${gencodeProtFa} output/$@.psl output/$@.work
	diff expected/$@.psl output/$@.psl
	pslCheck -verbose=0 output/$@.psl

testInterproPonAbeCat1BlatAlign: mkout
	rm -rf output/$@.psl output/$@.work
	${interproProteinTranscriptAlign} CAT1 ${ponAbeCat1Meta} ${ponAbeCat1RnaFa} ${ponAbeCat1ProtFa} output/$@.psl output/$@.work
	diff expected/$@.psl output/$@.psl
	pslCheck -verbose=0 output/$@.psl

interproTests:  interproMapTests interproAnnotsToDecoratorsTests

interproMapTests: testInterproAnnotsMap testInterproPonAbeCat1AnnotsMap

testInterproAnnotsMap: mkout
	${interproAnnotsMap} ${gencodePcPsl} \
	    --annot2TransPsl=output/$@.annotsTrans.psl \
	    ${gencodeInterproTsv} expected/testInterproBlatAlign.psl \
            output/$@.psl output/$@.ref.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.ref.tsv output/$@.ref.tsv
	diff expected/$@.annotsTrans.psl output/$@.annotsTrans.psl
	pslCheck -verbose=0 output/$@.psl

testInterproPonAbeCat1AnnotsMap: mkout
	${interproAnnotsMap} ${ponAbeCat1Psl} \
	    ${ponAbeCat1InterproTsv} expected/testInterproPonAbeCat1BlatAlign.psl \
            output/$@.psl output/$@.ref.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.ref.tsv output/$@.ref.tsv
	pslCheck -verbose=0 output/$@.psl

interproAnnotsToDecoratorsTests: testInterproAnnotsToDeco testInterproPonAbeCat1AnnotsToDeco

testInterproAnnotsToDeco: mkout
	${interproAnnotsToDecorators} --nproc=${nproc} ${gencodeInterproTsv} \
	    expected/testInterproAnnotsMap.psl expected/testInterproAnnotsMap.ref.tsv \
            output/$@.bed --featTypesTsv=output/$@.types.tsv
	diff expected/$@.bed output/$@.bed
	diff expected/$@.types.tsv output/$@.types.tsv
	bedToBigBed -type=bed12+ -as=${interproDecoAs} -tab output/$@.bed ${hg38ChromSizes} output/$@.bb

testInterproPonAbeCat1AnnotsToDeco: mkout
	${interproAnnotsToDecorators} --nproc=${nproc} ${ponAbeCat1InterproTsv} \
	    expected/testInterproPonAbeCat1AnnotsMap.psl expected/testInterproPonAbeCat1AnnotsMap.ref.tsv \
            output/$@.bed --featTypesTsv=output/$@.types.tsv
	diff expected/$@.bed output/$@.bed
	diff expected/$@.types.tsv output/$@.types.tsv
	bedToBigBed -type=bed12+ -as=${interproDecoAs} -tab output/$@.bed ${ponAbeChromSizes} output/$@.bb

##
# xspecies analysze tests
## 
xspeciesAnalyzeTests: testXsAnalyzePonAbe

testXsAnalyzePonAbe: mkout
	${xspeciesAnalyze} \
	    ${swissprotAnnots} \
	    expected/testXsPonAbeUniprotAnnotsMap.psl expected/testXsPonAbeUniprotAnnotsMap.ref.tsv \
	    CAT1 ${ponAbeCat1Psl} \
	    ${ponAbeCat1InterproTsv} \
	    expected/testInterproPonAbeCat1AnnotsMap.psl expected/testInterproPonAbeCat1AnnotsMap.ref.tsv \
	    ${uniprotInterproAnnotAssocTsv} output/$@.report.tsv


##
# Gene tracks for hub
##
bigGenePredAs = output/bigGenePred.as
gencodeBGpIn = output/gencode.v43.bpgin
gencodeGeneNames = output/gencode.v43.genenames
${gencodeBb}: ${gencodeGp} ${bigGenePredAs} mkout
	tawk '{print $$4,$$3,$$1}' ${gencodeMeta} > ${gencodeGeneNames}
	genePredToBigGenePred -geneNames=${gencodeGeneNames} $< /dev/stdout | sort -k1,1 -k2,2n > ${gencodeBGpIn}
	bedToBigBed -type=bed12 -tab -extraIndex=name,name2,geneName,geneName2 -as=${bigGenePredAs} ${gencodeBGpIn} ${hg38ChromSizes} $@

${bigGenePredAs}: mkout
	wget -O $@ -nv https://genome.ucsc.edu/goldenpath/help/examples/bigGenePred.as

${ponAbeCat1Bb}: ${ponAbeCat1BigCat} ${bigCatAs} mkout
	bedToBigBed -tab -type=bed12+ -as=${bigCatAs} -extraIndex=name,name2,geneName -sizesIs2Bit $< ${ponAbeTwobit} $@.tmp
	mv -f $@.tmp $@

${symSynCat1Bb}: ${symSynCat1BigCat} ${bigCatAs} mkout
	bedToBigBed -tab -type=bed12+ -as=${bigCatAs}  -extraIndex=name,name2,geneName -sizesIs2Bit $< ${symSynTwobit} $@.tmp
	mv -f $@.tmp $@

##
# hub
##
hub: ${gencodeBrowserDirDone} ${testHubFilesOutput} ${gencodeBb} ${ponAbeCat1Bb} ${symSynCat1Bb}
	hubCheck $$(pwd)/output/hub.txt

${gencodeBrowserDirDone}: ${gencodeGp}
	@mkdir -p ${gencodeBrowserDir}
	genePredToBed ${gencodeGp} /dev/stdout | \
	   bedToHtmlDir -hub ${testHubUrl} -title 'GENCODE Deco Test' hg38 /dev/stdin ${gencodeBrowserDir}
	touch $@

output/%: hub-files/%
	@mkdir -p $(dir $@)
	cp -f $< $@

mkout:
	@mkdir -p output


clean:
	rm -rf output


savebak:
	cd .. && ${MAKE} savebak


lint:
	cd .. && ${MAKE} lint
