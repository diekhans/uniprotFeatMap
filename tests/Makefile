proteinTranscriptAlign = ../bin/proteinTranscriptAlign 
uniprotGencodeSelect = ../bin/uniprotGencodeSelect

all:
	@echo "Note: targets are:"
	@echo "      make test"
	@echo "      make fulltest"
	@echo "      make aligntest"
	@echo "      make clean"


test: uniprotGencodeSelectTests
	@echo "to test blast alignments use: make fulltest" >&2

fulltest: test testRunBlastJobs
aligntest: testRunBlastJobs

uniprotGencodeSelectTests: testUniprotGencodeSelectBasic

testUniprotGencodeSelectBasic: mkout
	${uniprotGencodeSelect} input/gencode.v43.metadata.tsv input/swissprot.9606.tab expected/testRunBlastJobs.psl output/$@.psl output/$@.problems.tsv
	diff expected/$@.psl output/$@.psl
	diff expected/$@.problems.tsv output/$@.problems.tsv


testRunBlastJobs: mkout
	rm -rf output/$@.psl output/$@.work
	${proteinTranscriptAlign} input/swissprot.9606.fa input/gencode.v43.pc_transcripts.fa output/$@.psl output/$@.work
	cut -f 14 output/$@.work/aligns/*.psl | sort -u > output/$@.raw.tids
	cut -f 14 output/$@.psl | sort -u > output/$@.filtered.tids
	diff expected/$@.raw.tids output/$@.raw.tids
	diff expected/$@.filtered.tids output/$@.filtered.tids

mkout:
	@mkdir -p output


clean:
	rm -rf output

