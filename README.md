## mapping flow

annotation -> uniprotCanonical -> transcript -> genome


## Commands:

1. uniprotProteinTranscriptAlign - Align protein sequences to transcript RNAs with BLAT or BLAST. Creates a protein to NA PSL alignments, with some basic filtering.
   * input: uniprotFa, transFa
   * output: protTransPsl
1. uniprotProteinTranscriptMap - Filter protein to PSL transcript alignments to pair them based on being the listed transcript in UniProt.  Project the primary alignments to other transcript isoforms using the genomic coordiates.  This proved more accurate than doing protein alignments to other isoforms
   * input: gencodeMetaTsv, uniprotMetaTsv, protTransPsl
   * output: protTransPairedPsl, problemLogTsvh
1. uniprotMapAnnots - Map Uniprot annotations to the genome via protein and transcript alignments.  The output will be a NA to NA PSL alignments of annotations of all annotation types that are mapped.  They can be filtered later when building decorators.
   * input: uniprotAnnotsTsv, protTransPairedPsl, transGenomePsl
   * output: annotGenomePsl, annotTransRefTsv
1. uniprotAnnotsToDecorators - Convert domain annotations alignments create by uniprotMapAnnots to a decorator BED file in uniprotDecoration.as format.  Possibly filtering the results.
   * input: uniprotAnnotsTsv, annotGenomePsl, annotTransRefTsv
   * output: annotDecoratorBedFile
