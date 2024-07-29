# Mapping UniProt annotations to transcripts

## mapping flow

annotation -> uniprotCanonical -> transcript -> genome

## Terminology
- CDS to transcript: A protein query seqeunce converted to CDS coordinates (multiplied by 3) in alignment to transcript sequence so spliced codons can be tracked.
- canonical and non-canonical transcripts: Canonical transcripts are those listed by UniProt as matching a UniProt protein,.
- non-canonical transcripts: Transcript were the CDS overlaps the canonical transcript in the genome and we map the UniProt
- uniprot acc and id - The acc is the accession of a transcript, less the version, id includes the version.
- disrupted feature (disruption) - A disruption is a UniProt feature that only partially mapped to a non-canonical transcripts.
- deleted feature (deletion) - A deletion is a UniProt feature that doesn't map to a non-canonical transcripts.
- annotMapId: unique id assigned to the mapping of a UniProt feature annotation to a transcript, which is described below.

## Mapping pipeline commands:

1. `uniprotProteinTranscriptAlign`: Align protein sequences to transcript RNAs with BLAT or BLAST. Creates a protein to NA PSL alignments, with some basic filtering.
   * input: uniprotFa, uniprotMetaTsv, transFa
   * output: protCanonTransPsl
1. `uniprotProteinTranscriptMap`: Filter protein to PSL transcript alignments to pair them based on being the listed transcript in UniProt.  Project the primary alignments to other transcript isoforms using the genomic coordinates.  This proved more accurate than doing protein alignments to other isoforms.  Output is in UniProt CDS to transcript alignments.
   * input: gencodeMetaTsv, uniprotMetaTsv, protCanonTransPsl
   * output: cdsTransPairedPsl, problemLogTsvh
1. `uniprotAnnotsMap`: Map Uniprot annotations to the genome via protein and transcript alignments.  The output will be a CDS to transcript (NA to NA) PSL alignments of annotations of all annotation types that are mapped.  They can be filtered later when building decorators.  This can also map to other assembly via way of transcript-transcript alignments.
   * input: uniprotAnnotsTsv, cdsTransPairedPsl, transGenomePsl
   * output: annotGenomePsl, annotTransRefTsv
1. `uniprotAnnotsToDecorators`: Convert domain annotations alignments create by uniprotAnnotsMap to a decorator BED file in uniprotDecoration.as format.  Possibly filtering the results.
   * input: uniprotAnnotsTsv, annotGenomePsl, annotTransRefTsv
   * output: annotDecoratorBedFile

## Data acquisition commands
- `transMapMRnas`: Use genomic projection alignments to produce a set of cross-species mRNA to mRNA alignments.
- `getCatXSpeciesPairs`: Obtain paired GENCODE source transcripts and CAT cross-species annotations.

# Data required 

* U  /hive/data/outside/otto/uniprot/

# Cross-species projection

Domains can be projected to other species by way of mRNA to mRNA alignments
from the source assembly mRNA to the target assembly mRNA (`srcTargetMRnaPsl`).
This alignments is built and filtered based on the gene annotations

Source mRNAs (srcTrans) and aligned to target mRNAs (targetTrans) using
the TransMap algorithm, giving mRNA/mRNA alignments.  These are then
filter to pick the best transcript ortholog.  

For CAT annotations produced from GENCODE, this matches the transcript
accessions.  Other gene sets will require different criteria.

These alignments of srcTrans, annotated by `uniprotProteinTranscriptMap`, to
targetTrans in the target genomes, along with alignments of the target mRNAs
to the target genome are required.  These two alignments are given as
additional input to `uniprotAnnotsMap`, to project UniProt features through to
the target transcripts on the target genome.


# Identification

Each feature mapping needs to be assigned a unique id for use in accessing
individual feature annotations mappings, called the *annotation mapping id*
(`annotMapId`) and it is not stable multiple builds.  The ``annotMapId
is used for mouse clicks and other record identification.  It is the name column in
a decorator BED.

The format of a `annotMapId` is `<canon_acc>|<annot_idx>|<map_idx>`, where

- `canon_acc`: UniProt canonical isoform id
- `annot_idx`: row index of the feature in the UniProt annotation tab file, relative to the first feature for the accession
- `map_idx`: index of assignment of this feature annotation to a transcript mapping.  The number does not identify the same transcript across multiple features.

A file `*.ref.tsv` file is created when annotations are mapped to using in mapping `annotMapId`s back to the UniProt feature and transcript.

Disruption decorators will have an addition disruption index appended to `annotMapId` to make them unique.


