SHELL = /bin/bash -beEu -o pipefail

# used nprocs= to change number of processes
ifeq (${nprocs},)
  nprocs = 8
endif

uniprotProteinTranscriptAlign = ../bin/uniprotProteinTranscriptAlign 
uniprotGencodeSelect = ../bin/uniprotGencodeSelect
uniprotMapAnnots = ../bin/uniprotMapAnnots
uniprotAnnotsToDecorators = ../bin/uniprotAnnotsToDecorators

chromSizes = /hive/data/genomes/hg38/chrom.sizes
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


test: uniprotGencodeSelectTests uniprotMapAnnotsTests uniprotAnnotsToDecoratorsTests
	@echo "Note to test alignments use: make fulltest" >&2

fulltest: test aligntest
aligntest: protTransAlignTests

uniprotGencodeSelectTests: testUniprotGencodeSelectBasic

testUniprotGencodeSelectBasic: mkout
	${uniprotGencodeSelect} input/gencode.v43.metadata.tsv input/swissprot.9606.tab expected/testProtTransBlastAlign.psl output/$@.psl output/$@.problems.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.problems.tsv output/$@.problems.tsv
	pslCheck -verbose=0 output/$@.psl

define checkAlignments
	diff expected/$@.psl output/$@.psl
	pslCheck -verbose=0 output/$@.psl
endef

protTransAlignTests: testProtTransBlastAlign testProtTransBlatAlign
testProtTransBlastAlign: mkout
	rm -rf output/$@.psl output/$@.work
	${uniprotProteinTranscriptAlign} --algo=blast input/swissprot.9606.fa input/gencode.v43.pc_transcripts.fa output/$@.psl output/$@.work
	$(call checkAlignments)

testProtTransBlatAlign: mkout
	rm -rf output/$@.psl output/$@.work
	${uniprotProteinTranscriptAlign} --algo=blat input/swissprot.9606.fa input/gencode.v43.pc_transcripts.fa output/$@.psl output/$@.work
	$(call checkAlignments)

uniprotMapAnnotsTests: testUniprotMapAnnots

# intermediates saves for checking
testUniprotMapAnnots: mkout
	${uniprotMapAnnots} --interPrefix=output/$@.inter. \
	    input/swissprot.9606.annots.tab \
	    expected/testUniprotGencodeSelectBasic.psl \
	    input/gencode.v43.psl \
            output/$@.psl output/$@.ref.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.ref.tsv output/$@.ref.tsv
	diff expected/$@.inter.annotTrans.psl output/$@.inter.annotTrans.psl
	pslCheck -verbose=0 output/$@.psl
	pslCheck -verbose=0 output/$@.inter.annotTrans.psl

uniprotAnnotsToDecoratorsTests: testUniprotAnnotsToDecorators

testUniprotAnnotsToDecorators: mkout output/gencode.v43.bb ${testHubFilesOutput}
	${uniprotAnnotsToDecorators} --nprocs=${nprocs} input/swissprot.9606.tab input/swissprot.9606.annots.tab \
	    expected/testUniprotMapAnnots.psl expected/testUniprotMapAnnots.ref.tsv \
            input/gencode.v43.psl output/$@.bed --featTypesTsv=output/$@.types.tsv
	bedToBigBed -type=bed12+ -as=${decoAs} -tab output/$@.bed ${chromSizes} output/$@.bb
	diff expected/$@.bed output/$@.bed
	diff expected/$@.types.tsv output/$@.types.tsv

# gene track can can be used to validate decorators
output/gencode.v43.bb: input/gencode.v43.gp
	genePredToBed $< /dev/stdout | sort -k1,1 -k2,2n >output/gencode.v43.bed
	bedToBigBed -type=bed12 -tab output/gencode.v43.bed ${chromSizes} $@

output/%: hub-files/%
	mkdir -p $(dir $@)
	cp -f $< $@

mkout:
	@mkdir -p output


clean:
	rm -rf output

