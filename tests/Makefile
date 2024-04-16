root = ..
include ${root}/defs.mk

# see bin/getTestData for updating test data

# used nprocs= to change number of processes
ifeq (${nprocs},)
  nprocs = 8
endif

chromSizes = input/hg38.sizes
decoAs = ../etc/uniprotDecoration.as

testHubFiles = empty.html genomes.txt hub.txt trackDb.txt
testHubFilesOutput = ${testHubFiles:%=output/%}

all:
	@echo "Note: targets are:"
	@echo "      make test"
	@echo "      make fulltest"
	@echo "      make aligntest"
	@echo "      make clean"
	@exit 1


test: uniprotProteinTranscriptMapTests uniprotMapAnnotsTests uniprotAnnotsToDecoratorsTests uniprotDecoratorsMergeTests
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
	${uniprotProteinTranscriptAlign} --algo=${1} input/${2}.9606.fa input/gencode.v43.pc_transcripts.fa output/$@.psl output/$@.work
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

# $call runUniprotProteinTranscriptMap,swissprot,SP)
define runUniprotProteinTranscriptMap
	${uniprotProteinTranscriptMap} --interPrefix=output/$@.inter. \
		input/gencode.v43.metadata.tsv input/gencode.v43.gp input/gencode.v43.psl \
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
	    input/gencode.v43.psl \
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
            input/gencode.v43.psl output/$@.bed --featTypesTsv=output/$@.types.tsv
	diff expected/$@.bed output/$@.bed
	diff expected/$@.types.tsv output/$@.types.tsv
endef

testUniprotAnnotsToDecoratorsSP: mkout output/gencode.v43.bb ${testHubFilesOutput}
	$(call runUniprotAnnotsToDecorators,swissprot,SP)

testUniprotAnnotsToDecoratorsTR: mkout output/gencode.v43.bb ${testHubFilesOutput}
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
output/gencode.v43.bb: input/gencode.v43.gp
	genePredToBed $< /dev/stdout | sort -k1,1 -k2,2n >output/gencode.v43.bed
	bedToBigBed -type=bed12 -tab output/gencode.v43.bed ${chromSizes} $@

output/%: hub-files/%
	@mkdir -p $(dir $@)
	cp -f $< $@

mkout:
	@mkdir -p output


clean:
	rm -rf output

