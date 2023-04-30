SHELL = /bin/bash -beEu -o pipefail

proteinTranscriptAlign = ../bin/proteinTranscriptAlign 
uniprotGencodeSelect = ../bin/uniprotGencodeSelect
uniprotIsoCanonicalAlign = ../bin/uniprotIsoCanonicalAlign
uniprotIsoCanonicalSelect = ../bin/uniprotIsoCanonicalSelect

all:
	@echo "Note: targets are:"
	@echo "      make test"
	@echo "      make fulltest"
	@echo "      make aligntest"
	@echo "      make clean"


test: uniprotGencodeSelectTests uniprotIsoCanonicalSelectTests
	@echo "Note to test alignments use: make fulltest" >&2

fulltest: test aligntest
aligntest: testProtTransAlign testProtCanonAlign

uniprotGencodeSelectTests: testUniprotGencodeSelectBasic

testUniprotGencodeSelectBasic: mkout
	${uniprotGencodeSelect} input/gencode.v43.metadata.tsv input/swissprot.9606.tab expected/testProtTransAlign.psl output/$@.psl output/$@.problems.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.problems.tsv output/$@.problems.tsv

uniprotIsoCanonicalSelectTests: testUniprotIsoCanonicalSelectBasic

testUniprotIsoCanonicalSelectBasic: mkout
	${uniprotIsoCanonicalSelect} input/swissprot.9606.tab expected/testProtCanonAlign.psl output/$@.psl output/$@.problems.tsv
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

testProtTransAlign: mkout
	rm -rf output/$@.psl output/$@.work
	${proteinTranscriptAlign} input/swissprot.9606.fa input/gencode.v43.pc_transcripts.fa output/$@.psl output/$@.work
	$(call checkAlignments)

testProtCanonAlign: mkout
	rm -rf output/$@.psl output/$@.work
	${uniprotIsoCanonicalAlign} input/swissprot.9606.fa output/$@.psl output/$@.work
	$(call checkAlignments)


mkout:
	@mkdir -p output


clean:
	rm -rf output

