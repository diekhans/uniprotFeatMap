SHELL = /bin/bash -beEu -o pipefail

# used nprocs= to change number of processes
ifeq (${nprocs},)
  nprocs = 8
endif

proteinTranscriptAlign = ../bin/proteinTranscriptAlign 
uniprotGencodeSelect = ../bin/uniprotGencodeSelect
uniprotIsoCanonicalAlign = ../bin/uniprotIsoCanonicalAlign
uniprotIsoCanonicalSelect = ../bin/uniprotIsoCanonicalSelect
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


test: uniprotGencodeSelectTests uniprotIsoCanonicalSelectTests uniprotMapAnnotsTests uniprotAnnotsToDecoratorsTests
	@echo "Note to test alignments use: make fulltest" >&2

fulltest: test aligntest
aligntest: protTransAlignTests protCanonAlignTest

uniprotGencodeSelectTests: testUniprotGencodeSelectBasic

testUniprotGencodeSelectBasic: mkout
	${uniprotGencodeSelect} input/gencode.v43.metadata.tsv input/swissprot.9606.tab expected/testProtTransBlastAlign.psl output/$@.psl output/$@.problems.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.problems.tsv output/$@.problems.tsv
	pslCheck -verbose=0 output/$@.psl

uniprotIsoCanonicalSelectTests: testUniprotIsoCanonicalSelectBasic

testUniprotIsoCanonicalSelectBasic: mkout
	${uniprotIsoCanonicalSelect} input/swissprot.9606.tab expected/testProtCanonBlastAlign.psl output/$@.psl output/$@.problems.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.problems.tsv output/$@.problems.tsv

define checkAlignments
	pslStats -tsv -queries=input/swissprot.9606.sizes output/$@.psl output/$@.stats.tsv
	diff expected/$@.stats.tsv output/$@.stats.tsv
	diff expected/$@.psl output/$@.psl
	pslCheck -verbose=0 output/$@.psl
endef

protTransAlignTests: testProtTransBlastAlign testProtTransBlatAlign
testProtTransBlastAlign: mkout
	rm -rf output/$@.psl output/$@.work
	${proteinTranscriptAlign} --algo=blast input/swissprot.9606.fa input/gencode.v43.pc_transcripts.fa output/$@.psl output/$@.work
	$(call checkAlignments)

testProtTransBlatAlign: mkout
	rm -rf output/$@.psl output/$@.work
	${proteinTranscriptAlign} --algo=blat input/swissprot.9606.fa input/gencode.v43.pc_transcripts.fa output/$@.psl output/$@.work
	$(call checkAlignments)

protCanonAlignTest: testProtCanonBlastAlign testProtCanonBlatAlign
testProtCanonBlastAlign: mkout
	rm -rf output/$@.psl output/$@.work
	${uniprotIsoCanonicalAlign} --algo=blast input/swissprot.9606.fa output/$@.psl output/$@.work
	$(call checkAlignments)
testProtCanonBlatAlign: mkout
	rm -rf output/$@.psl output/$@.work
	${uniprotIsoCanonicalAlign} --algo=blat input/swissprot.9606.fa output/$@.psl output/$@.work
	$(call checkAlignments)


uniprotMapAnnotsTests: testUniprotMapAnnots

# intermediates saves for checking
testUniprotMapAnnots: mkout
	${uniprotMapAnnots} --interPrefix=output/$@.inter. \
	    input/swissprot.9606.annots.tab \
	    expected/testProtCanonBlatAlign.psl \
	    expected/testUniprotGencodeSelectBasic.psl \
	    input/gencode.v43.psl \
            output/$@.psl output/$@.ref.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.ref.tsv output/$@.ref.tsv
	diff expected/$@.inter.annotCanon.psl output/$@.inter.annotCanon.psl
	diff expected/$@.inter.annotIso.psl output/$@.inter.annotIso.psl
	diff expected/$@.inter.annotTrans.psl output/$@.inter.annotTrans.psl
	pslCheck -verbose=0 output/$@.psl
	pslCheck -verbose=0 output/$@.inter.annotCanon.psl
	pslCheck -verbose=0 output/$@.inter.annotIso.psl
	pslCheck -verbose=0 output/$@.inter.annotTrans.psl

uniprotAnnotsToDecoratorsTests: testUniprotAnnotsToDecorators

testUniprotAnnotsToDecorators: mkout output/gencode.v43.bb ${testHubFilesOutput}
	${uniprotAnnotsToDecorators} --nprocs=${nprocs} input/gencode.v43.psl input/swissprot.9606.annots.tab \
	    expected/testUniprotMapAnnots.psl expected/testUniprotMapAnnots.ref.tsv \
            output/$@.bed
	diff expected/$@.bed output/$@.bed
	bedToBigBed -type=bed12+ -as=${decoAs} -tab output/$@.bed ${chromSizes} output/$@.bb

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

