# Mapping UniProt and InterProScan annotations to transcripts

This modules implements mapping of protein annotations to transcripts and
generating BED track decorators.

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
- raw.psl - intermediate alignments that have not yet been fully filtered

## UniProt Mapping pipeline commands:

1. `uniprotProteinTranscriptAlign`: Align UniProt canonical protein sequences to transcript RNAs with BLAT or BLAST. Creates protein to NA PSL alignments for the transcripts associated with the UniProt canonical protein.
   * input: uniprotFa, uniprotMetaTsv, transFa
   * output: protCanonTransPsl
1. `uniprotProteinTranscriptMap`: Project the canonical alignments to other transcript isoforms using the genomic coordinates.  This proved more accurate than doing protein alignments to other isoforms.  Output is in the combination of the canonical and projected transcript alignments, with  the protein convert to CDS coordinates.
   * input: gencodeMetaTsv, gencodeGp, gencodePsl, uniprotMetaTsv, protCanonTransPsl
   * output: protCdsTransPairedPsl, problemLogTsv
1. `uniprotAnnotsMap`: Map Uniprot annotations to the genome via protein and transcript alignments.  The output will be a CDS to transcript (NA to NA) PSL alignments of annotations of all annotation types that are mapped.  They can be filtered later when building decorators.  This can also map to other assembly via way of transcript-transcript alignments.
   * input: uniprotAnnotsTsv, protCdsTransPairedPsl, transGenomePsl
   * output: annotGenomePsl, annotTransRefTsv
1. `uniprotAnnotsToDecorators`: Convert domain annotations alignments create by uniprotAnnotsMap to a decorator BED file in uniprotDecoration.as format.  Possibly filtering the results.
   * input: uniprotAnnotsTsv, annotGenomePsl, annotTransRefTsv
   * output: annotDecoratorBedFile

## Data acquisition commands
- `transMapMRnas`: Use genomic projection alignments to produce a set of cross-species mRNA to mRNA alignments.
- `getCatXSpeciesPairs`: Obtain paired GENCODE source transcripts and CAT cross-species annotations.

# Data required 

## UniProt data
* /hive/data/outside/otto/uniprot/

## Gene set specification

Due to different in gene set metadata, there is a requirement for
code for each supported gene set.

The data required for each gene set:
- metadata - format varies
- protein FASTA - protein sequences, metadata maps protein to transcript id.
  They can be obtain with genePredToProt or other methods.
- mRNA FASTA - CDS in upper case, UTR in lower case.
  They can be obtained with getRnaPred -cdsUpper.
- annotations to genome PSL alignments - which can be generated with genePredToPsl for annotations that don't have alignments.
- genePred format annotations - must exactly match the PSL.
  This is not needed for InterProScan.

The inputs file maybe compressed.

Different gene sets that are currently supported:

- GENCODE - GENCODE using the UCSC import formats
  - metadata in TSV format from UCSC import (not the table), which includes tags

- CAT1 - CAT used in CHM13 and initial primate references
  This is problematic in the fact that the name used the GENCODE transcript name, which is not stable and will not be unique when CAT maps paralogs.


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


# Uniprot Mapped Annotation Identification

Each feature mapping needs to be assigned a unique id for use in accessing
individual feature annotations mappings, called the *annotation mapping id*
(`annotMapId`) and it is not stable multiple builds.  The ``annotMapId
is used for mouse clicks and other record identification.  It is the name column in
a decorator BED.

The format of a `annotMapId` is `<canon_acc>|<annot_idx>|<map_idx>`, where

- `canon_acc`: UniProt canonical protein isoform id
- `annot_idx`: row index of the feature in the UniProt annotation tab file, relative to the first feature for the accession
- `map_idx`: index of assignment of this feature annotation to a transcript mapping, uniquely identifying the mapping.  Deleted features that don't map still have map_idx assigned to them.

A file `*.ref.tsv` file is created when annotations are mapped to using in mapping `annotMapId`s back to the UniProt feature and transcript.

Disruption decorators will have an addition disruption index appended to `annotMapId` to make them unique.


# InterProScan

1. `interproProteinTranscriptAlign`: Align InterProScan annotated protein sequences to transcript RNAs with BLAT or BLAST. Creates protein to NA PSL alignments for the transcript paired in the geneset.
   * input: gencodeMetaTsv, proteinFa, transFa
   * output: protTransPsl

## InterProScan Mapped Annotation Identification

The format of a `annotMapId` matches UniProt annotMapId, only using the
protein id used in the InterProScan run instead of the `canon_acc`.


A file `*.ref.tsv` file is created when annotations are mapped to using in mapping `annotMapId`s back to the UniProt feature and transcript.



