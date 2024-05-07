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


# Identification

Each feature mapping needs to be assigned a unique id for use in accessing
individual feature annotations mappings, called the *annotation mapping id*
(`annotMapId`) and it is not stable across multiple builds.  The ``annotMapId
is used for mouse clicks and other record identification.  It is the name column in
a decorator BED.

The format of a `annotMapId` is `<canon_acc>|<annot_idx>|<map_idx>`, where

- `canon_acc` - UniProt canonical isoform id
- `annot_idx` - row index of the feature in the UniProt annotation tab file, relative to the first feature for the accession
- `map_idx` - index of assignment of this feature annotation to a transcript mapping.  The number does not identify the same transcript across multiple feeatures.

A file `*.ref.tsv` file is created when annotations are mapped to using in mapping `annotMapId`s back to the UniProt feature and transcript.




