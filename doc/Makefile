
code = ./bin/mkFlowGraph ./bin/graphlib.py

parts = prot_trans_map uniprot_annot xspecies_trans_map xspecies_annot interpro_annot xspecies_compare
all_pdfs = ${parts:%=output/%.pdf}
flow_pdf = output/flow.pdf
# page = --page=8.5,11

all: ${flow_pdf}

${flow_pdf}: ${all_pdfs}
	@mkdir -p $(dir $@)
	pdfunite ${all_pdfs} $@.tmp.pdf
	mv -f $@.tmp.pdf $@

output/%.pdf: ${code}
	@mkdir -p $(dir $@)
	./bin/mkFlowGraph ${page} $* $@.tmp.pdf
	mv -f $@.tmp.pdf $@



clean:
	rm -rf output
