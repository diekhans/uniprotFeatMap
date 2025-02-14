root = ..
include ${root}/defs.mk

# see bin/getTestData for updating test data

# used nproc= to change number of processes
ifeq (${nproc},)
  nproc = 8
endif

chromSizes = input/hg38.sizes

testHubFiles = empty.html genomes.txt hub.txt hg38.trackDb.txt GCA_028885655.2.trackDb.txt
testHubFilesOutput = ${testHubFiles:%=output/%}

testHubUrl = https://hgwdev.gi.ucsc.edu/~markd/kznf/uniprotFeatMap/tests/output/hub.txt

gencodeGp = input/gencode.v43.pc.gp
gencodeMeta = input/gencode.v43.metadata.tsv
gencodeBb = output/gencode.v43.bb
gencodeBrowserDir = output/gencode-dir
gencodeBrowserDirDone = ${gencodeBrowserDir}/done

ponAbeCatBb = output/GCA_028885655.2.cat.bb
ponAbeTwobit = https://hgdownload.soe.ucsc.edu/hubs/GCA/028/885/655/GCA_028885655.2/GCA_028885655.2.2bit



all:
	@echo "Note: targets are:"
	@echo "      make test"
	@echo "      make fulltest"
	@echo "      make aligntest"
	@echo "      make hub"
	@echo "      make clean"
	@exit 1


test: uniprotTests interproTests hub
	@echo "==========================================" >&2
	@echo "Note to test alignments use: make fulltest" >&2
	@echo "==========================================" >&2

fulltest: aligntest test

# alignment tests must be done serially due to parasol losing small jobs on multiple
# batches with short jobs
aligntest:
	${MAKE} protTransAlignTests
	${MAKE} interproAlignTests

####
uniprotTests: uniprotProteinTranscriptMapTests uniprotAnnotsMapTests uniprotAnnotsToDecoratorsTests \
	uniprotDecoratorsMergeTests xspeciesTests uniprotInfoTests 


# serial due to parasol bug
protTransAlignTests:
	${MAKE} testProtTransSPBlastAlign
	${MAKE} testProtTransSPBlatAlign
	${MAKE} testProtTransTRBlastAlign
	${MAKE} testProtTransTRBlatAlign

# $(call runUpAlignments,blat,swissprot)
define runUpAlignments
	rm -rf output/$@.psl output/$@.work
	${uniprotProteinTranscriptAlign} --algo=${1} input/${2}.9606.tab input/${2}.9606.fa input/gencode.v43.pc.fa output/$@.psl output/$@.work
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
		GENCODE ${gencodeMeta} input/gencode.v43.pc.psl ${gencodeGp} \
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
	    input/gencode.v43.pc.psl \
	    input/${1}.9606.annots.tab expected/testUniprotProteinTranscriptMap${2}.psl \
            output/$@.psl output/$@.ref.tsv output/$@.problems.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.ref.tsv output/$@.ref.tsv
	diff expected/$@.problems.tsv output/$@.problems.tsv
	diff expected/$@.inter.annotTrans.psl output/$@.inter.annotTrans.psl
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
	${uniprotAnnotsToDecorators} --nproc=${nproc} input/gencode.v43.pc.psl \
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

# error handling, use to not terminate on fail to open files
testUniprotAnnotsToDecoratorsError: mkout ${gencodeBb}
	! ${uniprotAnnotsToDecorators} --nproc=${nproc} input/gencode.v43.pc.psl \
	    input/swissprot.9606.tab input/swissport.9606.annots.tab \
	    expected/testUniprotAnnotsMapSP.psl expected/does-not-exist.tsv \
            output/$@.bed >& output/$@.err || { echo "Error: command should fail"; exit 1; }
	grep -q 'Exception: Pool initialization failed' output/$@.err 

###
uniprotDecoratorsMergeTests: testUniprotDecoratorsMerge

testUniprotDecoratorsMerge: mkout
	${uniprotDecoratorsMerge} --outBed=output/$@.bed expected/testUniprotAnnotsToDecoratorsSP.bed expected/testUniprotAnnotsToDecoratorsTR.bed
	bedToBigBed -type=bed12+ -as=${uniprotDecoAs} -tab output/$@.bed ${chromSizes} output/$@.bb
	diff expected/$@.bed output/$@.bed

###
# cross-species mapping
xspeciesTests: testXSpeciesRnaMap testXSpeciesGencodeCatFilter testXSpeciesUniprotAnnotsMap \
	testXSpeciesAnnotsToDeco

testXSpeciesRnaMap: mkout
	${xspeciesTrans2TransMap} input/gencode.v43.pc.psl input/gencode.v43.pc.fa \
		input/GCA_028885655.2.chains input/GCA_028885655.2.cat.psl input/GCA_028885655.2.cat.fa output/$@.psl
	diff expected/$@.psl output/$@.psl

testXSpeciesGencodeCatFilter: mkout ${ponAbeCatBb}
	${xspeciesGencode2CatFilter} expected/testXSpeciesRnaMap.psl ${ponAbeCatBb} output/$@.psl
	diff expected/$@.psl output/$@.psl

testXSpeciesUniprotAnnotsMap: mkout
	${uniprotAnnotsMap} --interPrefix=output/$@.inter. \
	    --xspeciesTransPsl=expected/testXSpeciesGencodeCatFilter.psl \
	    input/GCA_028885655.2.cat.psl \
	    input/swissprot.9606.annots.tab expected/testUniprotProteinTranscriptMapSP.psl \
            output/$@.psl output/$@.ref.tsv output/$@.problems.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.ref.tsv output/$@.ref.tsv
	diff expected/$@.problems.tsv output/$@.problems.tsv
	diff expected/$@.inter.annotTrans.psl output/$@.inter.annotTrans.psl
	pslCheck -verbose=0 output/$@.psl
	pslCheck -verbose=0 output/$@.inter.annotTrans.psl

testXSpeciesAnnotsToDeco: mkout
	${uniprotAnnotsToDecorators} --nproc=${nproc} input/GCA_028885655.2.cat.psl \
	    input/swissprot.9606.tab input/swissprot.9606.annots.tab \
	    expected/testXSpeciesUniprotAnnotsMap.psl expected/testXSpeciesUniprotAnnotsMap.ref.tsv \
            output/$@.bed --featTypesTsv=output/$@.types.tsv
	diff expected/$@.bed output/$@.bed
	diff expected/$@.types.tsv output/$@.types.tsv
	bedToBigBed -type=bed12+ -as=${uniprotDecoAs} -sizesIs2Bit -tab expected/$@.bed ${ponAbeTwobit} output/$@.bb

###
uniprotInfoTests: testUniprotInfoEnsTrans testUniprotInfoGeneName

testUniprotInfoEnsTrans: mkout
	${uniprotInfo} ensemblTransId input/swissprot.9606.tab output/$@.tsv
	diff expected/$@.tsv output/$@.tsv

testUniprotInfoGeneName:
	${uniprotInfo} geneName input/swissprot.9606.tab output/$@.tsv
	diff expected/$@.tsv output/$@.tsv

####
# $(call runUpAlignments,blat)
define runIpAlignments
	rm -rf output/$@.psl output/$@.work
	${interproProteinTranscriptAlign} --algo=${1} GENCODE ${gencodeMeta} input/gencode.v43.pc.fa input/gencode.v43.proteins.fa output/$@.psl output/$@.work
	diff expected/$@.psl output/$@.psl
	pslCheck -verbose=0 output/$@.psl
endef

# serial due to parasol bug
interproAlignTests:
	${MAKE} testInterproBlatAlign

testInterproBlatAlign: mkout
	$(call runIpAlignments,blat)

interproTests: testInterproAnnotsMap testInterproAnnotsToDeco

testInterproAnnotsMap: mkout
	${interproAnnotsMap} \
	    input/gencode.v43.pc.psl \
	    input/gencode.v43.interproscan.tsv expected/testInterproBlatAlign.psl \
            output/$@.psl output/$@.ref.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.ref.tsv output/$@.ref.tsv
	pslCheck -verbose=0 output/$@.psl

testInterproAnnotsToDeco: mkout
	${interproAnnotsToDecorators} --nproc=${nproc} input/gencode.v43.interproscan.tsv \
	    expected/testInterproAnnotsMap.psl expected/testInterproAnnotsMap.ref.tsv \
            output/$@.bed --featTypesTsv=output/$@.types.tsv
	diff expected/$@.bed output/$@.bed
	diff expected/$@.types.tsv output/$@.types.tsv
	bedToBigBed -type=bed12+ -as=${interproDecoAs} -tab output/$@.bed ${chromSizes} output/$@.bb

##
# Gene tracks for hub
##
bigGenePredAs = output/bigGenePred.as
gencodeBGpIn = output/gencode.v43.bpgin
gencodeGeneNames = output/gencode.v43.genenames
${gencodeBb}: ${gencodeGp} ${bigGenePredAs} mkout
	tawk '{print $$4,$$3,$$1}' ${gencodeMeta} > ${gencodeGeneNames}
	genePredToBigGenePred -geneNames=${gencodeGeneNames} $< /dev/stdout | sort -k1,1 -k2,2n > ${gencodeBGpIn}
	bedToBigBed -type=bed12 -tab -as=${bigGenePredAs} ${gencodeBGpIn} ${chromSizes} $@

${bigGenePredAs}: mkout
	wget -O $@ -nv https://genome.ucsc.edu/goldenpath/help/examples/bigGenePred.as

${ponAbeCatBb}: input/GCA_028885655.2.cat.bigCat input/bigCat.as mkout
	bedToBigBed -tab -type=bed12+ -as=input/bigCat.as -extraIndex=name -sizesIs2Bit $< ${ponAbeTwobit} $@.tmp
	mv -f $@.tmp $@

##
# hub
##
hub: ${gencodeBrowserDirDone} ${testHubFilesOutput} ${gencodeBb} ${ponAbeCatBb}
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
