proteinTranscriptAlign = ../bin/proteinTranscriptAlign 


all:
	@echo "Note: targets are:"
	@echo "      make test"
	@echo "      make fulltest"
	@echo "      make clean"


test:

fulltest: test testRunBlastJobs


testRunBlastJobs: mkout
	rm -rf output/$@.psl output/$@.work
	${proteinTranscriptAlign} input/swissprot.9606.fa input/gencode.v43.pc_transcripts.fa output/$@.psl output/$@.work
	diff expected/$@.psl output/$@.psl

mkout:
	@mkdir -p output


clean:
	rm -rf output

