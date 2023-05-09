SHELL = /bin/bash -beEu -o pipefail

proteinTranscriptAlign = ../bin/proteinTranscriptAlign 
uniprotGencodeSelect = ../bin/uniprotGencodeSelect
uniprotIsoCanonicalAlign = ../bin/uniprotIsoCanonicalAlign
uniprotIsoCanonicalSelect = ../bin/uniprotIsoCanonicalSelect
uniprotMapAnnots = ../bin/uniprotMapAnnots

all:
	@echo "Note: targets are:"
	@echo "      make test"
	@echo "      make fulltest"
	@echo "      make aligntest"
	@echo "      make clean"


test: uniprotGencodeSelectTests uniprotIsoCanonicalSelectTests uniprotMapAnnotsTests
	@echo "Note to test alignments use: make fulltest" >&2

fulltest: test aligntest
aligntest: protTransAlignTests protCanonAlignTest

uniprotGencodeSelectTests: testUniprotGencodeSelectBasic

testUniprotGencodeSelectBasic: mkout
	${uniprotGencodeSelect} input/gencode.v43.metadata.tsv input/swissprot.9606.tab expected/testProtTransBlastAlign.psl output/$@.psl output/$@.problems.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.problems.tsv output/$@.problems.tsv

uniprotIsoCanonicalSelectTests: testUniprotIsoCanonicalSelectBasic

testUniprotIsoCanonicalSelectBasic: mkout
	${uniprotIsoCanonicalSelect} input/swissprot.9606.tab expected/testProtCanonBlastAlign.psl output/$@.psl output/$@.problems.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.problems.tsv output/$@.problems.tsv

define checkAlignments
	pslStats -tsv -queries=input/swissprot.9606.sizes output/$@.psl output/$@.stats.tsv
	cut -f 10,14 output/$@.work/aligns/*.psl | sort -u -k2,2 -k1,1 > output/$@.raw.tids
	cut -f 10,14 output/$@.psl | sort -u -k2,2 -k1,1 > output/$@.filtered.tids
	diff expected/$@.raw.tids output/$@.raw.tids
	diff expected/$@.filtered.tids output/$@.filtered.tids
	diff expected/$@.stats.tsv output/$@.stats.tsv
	diff expected/$@.psl output/$@.psl
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

testUniprotMapAnnots: mkout
	${uniprotMapAnnots} input/swissprot.9606.annots.tab \
	    expected/testProtCanonBlatAlign.psl \
	    expected/testUniprotGencodeSelectBasic.psl \
	    input/gencode.v43.psl \
            output/$@.psl
	diff expected/$@.psl output/$@.psl

mkout:
	@mkdir -p output


clean:
	rm -rf output

