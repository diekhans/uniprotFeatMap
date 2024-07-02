root = ..
include ${root}/defs.mk

# see bin/getTestData for updating test data

# used nprocs= to change number of processes
ifeq (${nprocs},)
  nprocs = 8
endif

chromSizes = input/hg38.sizes

testHubFiles = empty.html genomes.txt hub.txt trackDb.txt
testHubFilesOutput = ${testHubFiles:%=output/%}

testHubUrl = https://hgwdev.gi.ucsc.edu/~markd/kznf/domainDecorators/tests/output/hub.txt

gencodeBed = output/gencode.v43.bed
gencodeBb = output/gencode.v43.bb
gencodeBrowserDir = output/gencode-dir
gencodeBrowserDirDone = ${gencodeBrowserDir}/done

all:
	@echo "Note: targets are:"
	@echo "      make test"
	@echo "      make fulltest"
	@echo "      make aligntest"
	@echo "      make buildHub"
	@echo "      make clean"
	@exit 0


test: uniprotProteinTranscriptMapTests uniprotMapAnnotsTests uniprotAnnotsToDecoratorsTests \
	uniprotDecoratorsMergeTests xspeciesTests, buildHub
	@echo "==========================================" >&2
	@echo "Note to test alignments use: make fulltest" >&2
	@echo "==========================================" >&2

fulltest: aligntest test
aligntest: protTransAlignTests

###
protTransAlignTests: testProtTransSPBlastAlign testProtTransSPBlatAlign \
	testProtTransTRBlastAlign testProtTransTRBlatAlign \

# $(call runAlignments,blat,swissprot)
define runAlignments
	rm -rf output/$@.psl output/$@.work
	${uniprotProteinTranscriptAlign} --algo=${1} input/${2}.9606.fa input/gencode.v43.pc.fa output/$@.psl output/$@.work
	diff expected/$@.psl output/$@.psl
	pslCheck -verbose=0 output/$@.psl
endef

testProtTransSPBlastAlign: mkout
	$(call runAlignments,blast,swissprot)

testProtTransSPBlatAlign: mkout
	$(call runAlignments,blat,swissprot)

testProtTransTRBlastAlign: mkout
	$(call runAlignments,blast,trembl)

testProtTransTRBlatAlign: mkout
	$(call runAlignments,blat,trembl)

####
uniprotProteinTranscriptMapTests: testUniprotProteinTranscriptMapSP testUniprotProteinTranscriptMapTR

# turn on for debugging
# uniprotProteinTranscriptMapInterPrefix = --interPrefix=output/$@.inter.
uniprotProteinTranscriptMapInterPrefix =

# $call runUniprotProteinTranscriptMap,swissprot,SP)
define runUniprotProteinTranscriptMap
	${uniprotProteinTranscriptMap} ${uniprotProteinTranscriptMapInterPrefix} \
		input/gencode.v43.metadata.tsv input/gencode.v43.pc.gp input/gencode.v43.pc.psl \
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
uniprotMapAnnotsTests: testUniprotMapAnnotsSP testUniprotMapAnnotsTR
# intermediates saved for debugging

# $(call runMapAnnots,swissprot,SP)
define runMapAnnots
	${uniprotMapAnnots} --interPrefix=output/$@.inter. \
	    input/${1}.9606.annots.tab expected/testUniprotProteinTranscriptMap${2}.psl \
	    input/gencode.v43.pc.psl \
            output/$@.psl output/$@.ref.tsv output/$@.problems.tsv
	$(call checkAnnots)
	diff expected/$@.psl output/$@.psl
	diff expected/$@.ref.tsv output/$@.ref.tsv
	diff expected/$@.problems.tsv output/$@.problems.tsv
	diff expected/$@.inter.annotTrans.psl output/$@.inter.annotTrans.psl
	pslCheck -verbose=0 output/$@.psl
	pslCheck -verbose=0 output/$@.inter.annotTrans.psl
endef

testUniprotMapAnnotsSP: mkout
	$(call runMapAnnots,swissprot,SP)
testUniprotMapAnnotsTR: mkout
	$(call runMapAnnots,trembl,TR)

###
uniprotAnnotsToDecoratorsTests: testUniprotAnnotsToDecoratorsSP testUniprotAnnotsToDecoratorsTR \
	testUniprotAnnotsToDecoColorHelp

# $(call runUniprotAnnotsToDecorators,swissprot,SP)
define runUniprotAnnotsToDecorators
	${uniprotAnnotsToDecorators} --nprocs=${nprocs} input/$(1).9606.tab input/$(1).9606.annots.tab \
	    expected/testUniprotMapAnnots$(2).psl expected/testUniprotMapAnnots$(2).ref.tsv \
            input/gencode.v43.pc.psl output/$@.bed --featTypesTsv=output/$@.types.tsv
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

###
uniprotDecoratorsMergeTests: testUniprotDecoratorsMerge

testUniprotDecoratorsMerge: mkout
	${uniprotDecoratorsMerge} --outBed=output/$@.bed expected/testUniprotAnnotsToDecoratorsSP.bed expected/testUniprotAnnotsToDecoratorsTR.bed
	bedToBigBed -type=bed12+ -as=${decoAs} -tab output/$@.bed ${chromSizes} output/$@.bb
	diff expected/$@.bed output/$@.bed


# gene track can can be used to validate decorators
${gencodeBb}: input/gencode.v43.pc.gp
	genePredToBed $< /dev/stdout | sort -k1,1 -k2,2n > ${gencodeBed}
	bedToBigBed -type=bed12 -tab ${gencodeBed} ${chromSizes} $@


# cross-species mapping
cat_bigcat = output/GCA_028885655.2.cat.bb
ponAbeTwobit = https://hgdownload.soe.ucsc.edu/hubs/GCA/028/885/655/GCA_028885655.2/GCA_028885655.2.2bit

xspeciesTests: testXSpeciesRnaMap testXSpeciesGencodeCatFilter testXSpeciesUniprotMapAnnots

testXSpeciesRnaMap: mkout
	${xspeciesTrans2TransMap} input/gencode.v43.pc.psl input/gencode.v43.pc.fa \
		input/GCA_028885655.2.chains input/GCA_028885655.2.cat.psl input/GCA_028885655.2.cat.fa output/$@.psl
	diff expected/$@.psl output/$@.psl

testXSpeciesGencodeCatFilter: mkout ${cat_bigcat}
	${xspeciesGencode2CatFilter} expected/testXSpeciesRnaMap.psl ${cat_bigcat} output/$@.psl
	diff expected/$@.psl output/$@.psl

testXSpeciesUniprotMapAnnots: mkout
	${uniprotMapAnnots} --interPrefix=output/$@.inter. \
	    --xspeciesTransPsl=expected/testXSpeciesGencodeCatFilter.psl \
	    input/swissprot.9606.annots.tab expected/testUniprotProteinTranscriptMapSP.psl \
	    input/GCA_028885655.2.cat.psl \
            output/$@.psl output/$@.ref.tsv output/$@.problems.tsv
	$(call checkAnnots)
	diff expected/$@.psl output/$@.psl
	diff expected/$@.ref.tsv output/$@.ref.tsv
	diff expected/$@.problems.tsv output/$@.problems.tsv
	diff expected/$@.inter.annotTrans.psl output/$@.inter.annotTrans.psl
	pslCheck -verbose=0 output/$@.psl
	pslCheck -verbose=0 output/$@.inter.annotTrans.psl

${cat_bigcat}: input/GCA_028885655.2.cat.bigCat input/bigCat.as
	@mkdir -p $(dir $@)
	bedToBigBed -tab -type=bed12+ -as=input/bigCat.as -extraIndex=name -sizesIs2Bit $< ${ponAbeTwobit} $@.tmp
	mv -f $@.tmp $@


# hub
buildHub: ${gencodeBrowserDirDone} ${testHubFilesOutput}

${gencodeBrowserDirDone}: ${gencodeBb}
	@mkdir -p ${gencodeBrowserDir}
	bedToHtmlDir -hub ${testHubUrl} -title 'GENCODE Deco Test' hg38 ${gencodeBed} ${gencodeBrowserDir}
	touch $@

output/%: hub-files/%
	@mkdir -p $(dir $@)
	cp -f $< $@

mkout:
	@mkdir -p output


clean:
	rm -rf output

