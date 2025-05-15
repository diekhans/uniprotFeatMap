
code = ./bin/mkFlowGraph ./bin/graphlib.py

all: flow xspecies

flow: output/mapping-flow.pdf
output/mapping-flow.pdf: ${code}
	@mkdir -p $(dir $@)
	./bin/mkFlowGraph $@.tmp.pdf
	mv -f $@.tmp.pdf $@

xspecies: output/xspecies-mapping-flow.pdf
output/xspecies-mapping-flow.pdf: ${code}
	@mkdir -p $(dir $@)
	./bin/mkFlowGraph --xspecies $@.tmp.pdf
	mv -f $@.tmp.pdf $@


clean:
	rm -rf output
