## mapping flow

annotation -> uniprotCanonical -> transcript -> genome

## Treminology
- CDS to transcript - protein query seqeunce converted to CDS coordinates (multiplied by 3) in alignment to trnascript sequecne so spliced codons can be tracked.
- canonical and non-canonical transcripts - canonical transcripts are those listed by UniProt as matching a particular protein.
- id and acc - acc is the accession of a transcript, less the version; id includes that version

## Commands:

1. uniprotProteinTranscriptAlign - Align protein sequences to transcript RNAs with BLAT or BLAST. Creates a protein to NA PSL alignments, with some basic filtering.
   * input: uniprotFa, transFa
   * output: protTransPsl
1. uniprotProteinTranscriptMap - Filter protein to PSL transcript alignments to pair them based on being the listed transcript in UniProt.  Project the primary alignments to other transcript isoforms using the genomic coordiates.  This proved more accurate than doing protein alignments to other isoforms.  Output is in UniProt CDS to transcri
   * input: gencodeMetaTsv, uniprotMetaTsv, protTransPsl
   * output: cdsTransPairedPsl, problemLogTsvh
1. uniprotMapAnnots - Map Uniprot annotations to the genome via protein and transcript alignments.  The output will be a CDS to transcript (NA to NA) PSL alignments of annotations of all annotation types that are mapped.  They can be filtered later when building decorators.
   * input: uniprotAnnotsTsv, cdsTransPairedPsl, transGenomePsl
   * output: annotGenomePsl, annotTransRefTsv
1. uniprotAnnotsToDecorators - Convert domain annotations alignments create by uniprotMapAnnots to a decorator BED file in uniprotDecoration.as format.  Possibly filtering the results.
   * input: uniprotAnnotsTsv, annotGenomePsl, annotTransRefTsv
   * output: annotDecoratorBedFile

# Data required 


gencode.v44.pc_transcripts.meta.tsv
gencode.v44.pc_transcripts.gp
gencode.v44.pc_transcripts.psl
gencode.v44.pc_transcripts.ids.tsv

swissprot.9606.fa.gz
swissprot.9606.refs.tab
swissprot.9606.annots.tab
swissprot.9606.tab

ln /hive/data/genomes/hg38/bed/gencodeV45/hgcImportPre/data/release_45/gencode.v45.pc_transcripts.fa.gz .
zcat gencode.v45.pc_transcripts.fa.gz | awk 'BEGIN{print "transcriptId"} /^>/{print gensub("^>([.A-Z0-9]+).*$", "\\1", 1)}' >  gencode.v45.pc_transcripts.ids.tsv

