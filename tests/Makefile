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


test: uniprotProteinTranscriptMapTests uniprotMapAnnotsTests uniprotAnnotsToDecoratorsTests
	@echo "Note to test alignments use: make fulltest" >&2

fulltest: test aligntest
aligntest: protTransAlignTests

###
protTransAlignTests: testProtTransSPBlastAlign testProtTransSPBlatAlign

testProtTransSPBlastAlign: mkout
	rm -rf output/$@.psl output/$@.work
	${uniprotProteinTranscriptAlign} --algo=blast input/swissprot.9606.fa input/gencode.v43.pc_transcripts.fa output/$@.psl output/$@.work
	$(call checkAlignments)

testProtTransSPBlatAlign: mkout
	rm -rf output/$@.psl output/$@.work
	${uniprotProteinTranscriptAlign} --algo=blat input/swissprot.9606.fa input/gencode.v43.pc_transcripts.fa output/$@.psl output/$@.work
	$(call checkAlignments)
####
uniprotProteinTranscriptMapTests: testUniprotProteinTranscriptMapSP

testUniprotProteinTranscriptMapSP: mkout
	${uniprotProteinTranscriptMap} --interPrefix=output/$@.inter. \
		input/gencode.v43.metadata.tsv input/gencode.v43.gp input/gencode.v43.psl \
		input/swissprot.9606.tab expected/testProtTransSPBlatAlign.psl output/$@.psl output/$@.problems.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.problems.tsv output/$@.problems.tsv
	pslCheck -verbose=0 output/$@.psl

define checkAlignments
	diff expected/$@.psl output/$@.psl
	pslCheck -verbose=0 output/$@.psl
endef

###
uniprotMapAnnotsTests: testUniprotMapAnnotsSP

# intermediates saved for debugging
testUniprotMapAnnotsSP: mkout
	${uniprotMapAnnots} --interPrefix=output/$@.inter. \
	    input/swissprot.9606.annots.tab \
	    expected/testUniprotProteinTranscriptMapSP.psl \
	    input/gencode.v43.psl \
            output/$@.psl output/$@.ref.tsv output/$@.problems.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.ref.tsv output/$@.ref.tsv
	diff expected/$@.problems.tsv output/$@.problems.tsv
	diff expected/$@.inter.annotTrans.psl output/$@.inter.annotTrans.psl
	pslCheck -verbose=0 output/$@.psl
	pslCheck -verbose=0 output/$@.inter.annotTrans.psl

###
uniprotAnnotsToDecoratorsTests: testUniprotAnnotsToDecoratorsSP \
	testUniprotAnnotsToDecoColorHelp

testUniprotAnnotsToDecoratorsSP: mkout output/gencode.v43.bb ${testHubFilesOutput}
	${uniprotAnnotsToDecorators} --nprocs=${nprocs} input/swissprot.9606.tab input/swissprot.9606.annots.tab \
	    expected/testUniprotMapAnnotsSP.psl expected/testUniprotMapAnnotsSP.ref.tsv \
            input/gencode.v43.psl output/$@.bed --featTypesTsv=output/$@.types.tsv
	bedToBigBed -type=bed12+ -as=${decoAs} -tab output/$@.bed ${chromSizes} output/$@.bb
	diff expected/$@.bed output/$@.bed
	diff expected/$@.types.tsv output/$@.types.tsv

testUniprotAnnotsToDecoColorHelp: mkout
	${uniprotAnnotsToDecorators} --help-colors >output/$@.out
	diff expected/$@.out output/$@.out

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

