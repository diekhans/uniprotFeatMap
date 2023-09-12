## commands


1. proteinTranscriptAlign - Align protein sequences to transcript RNAs with BLAT or BLAST. Creates a protein to NA PSL alignments, with some basic filtering.
   * input: protFa, transFa
   * output: protTransPsl
1. uniprotGencodeSelect - Filter protein to PSL transcript alignments to pair them based on gene.
   * input: gencodeMetaTsv, uniprotMetaTsv, protTransPsl
   * output: protTransPairedPsl, problemLogTsvh
1. uniprotIsoCanonicalAlign  - Produce mapping alignments of each UniProt isoform to the canonical isoform, including a 1-to-1 mapping of the canonical.  Result are protein to protein PSL alignments.
   * input: uniprotFa
   * output: isoCanonicalRawPsl
1. uniprotIsoCanonicalSelect - Filter alignments of UniProt proteins to keep only the alignments from to the canonical UniProt protein which has the annotations.
   * input: uniprotMetaTsv, isoCanonicalRawPsl
   * output: isoCanonicalPsl, problemLogTsv
1. uniprotMapAnnots - Map Uniprot annotations to the genome via protein and transcript alignments.  The output will be a NA to NA PSL alignments of annotations of all annotation types that are mapped.  They can be filtered later when building decorators.
   * input: uniprotAnnotsTsv, protTransPairedPsl, transGenomePsl
   * output: annotGenomePsl, annotTransRefTsv
1. uniprotAnnotsToDecorators - Convert domain annotations alignments create by uniprotMapAnnots to a decorator BED file in uniprotDecoration.as format.  Possibly filtering the results.
   * input: uniprotAnnotsTsv, annotGenomePsl, annotTransRefTsv
   * output: annotDecoratorBed
   
   
