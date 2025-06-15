"""
Annotation mappings data for analysis
"""
from collections import namedtuple, defaultdict
from pycbio.hgdata.psl import PslReader
from uniprotmap.metadata import annot2GenomeRefReader


class MappingError(Exception):
    pass


class AnnotMapping(namedtuple("AnnotMapping",
                              ("annotRef", "annotPsl", "annot"))):
    """A single mapping for an annotation to a genome. annotPsl is None if not mapped.
    The annot field depends on the type of annotation (UniProt or Interpro).
    """
    __slots__ = ()


class TransAnnotMappings(namedtuple("TransAnnotMappings",
                                    ("transcriptId", "chrom", "transPsl", "annotMappings",))):
    """Mappings for all annotations on a transcript."""
    __slots__ = ()


class AnnotMappingsTbl(list):
    """Table of TransAnnotMappings, all the mappings annotations to a genome"""
    # FIXME: still needed?
    def __init__(self, annot2GenomePslFile, annot2GenomeRefTsv):
        self.byTransId = defaultdict(list)
        for transAnnotMappings in transAnnotMappingReader(annot2GenomePslFile, annot2GenomeRefTsv):
            self._add(transAnnotMappings)
        self.byTransId.default_factory = None

    def _add(self, transAnnotMappings):
        self.byTransId[transAnnotMappings.transcriptId].append(transAnnotMappings)

    def findEntries(self, transId):
        return self.byTransId.get(transId, ())

    def getEntries(self, transId):
        entries = self.findEntries(transId)
        if len(entries) == 0:
            raise MappingError(f"transcript not found for '{transId}'")
        return entries

    def findEntry(self, transId, chrom):
        # handle multi-location entries
        entries = self.byTransId.get(transId, ())
        for entry in entries:
            if entry.chrom == chrom:
                return entry
        return None

    def getEntry(self, transId, chrom):
        entry = self.findEntry(transId, chrom)
        if entry is None:
            raise MappingError(f"transcript not found for '{transId}' on chrom '{chrom}'")
        return entry


###
# Reading annotation mappings
###
def _makeAnnotMapping(annot2GenomeRef, annot2GenomePsls, annot):
    if annot2GenomeRef.alignIdx is None:
        return AnnotMapping(annot2GenomeRef, None, annot)
    else:
        return AnnotMapping(annot2GenomeRef, annot2GenomePsls[annot2GenomeRef.alignIdx], annot)

def _makeTransAnnotMapping(annotMappings, transPslLookupFunc):
    annotRef0 = annotMappings[0].annotRef
    transPsl = None
    if transPslLookupFunc is not None:
        transPsl = transPslLookupFunc(annotRef0.transcriptId,
                                      annotRef0.transcriptPos.name)
    return TransAnnotMappings(annotRef0.transcriptId,
                              annotRef0.transcriptPos.name,
                              transPsl, tuple(annotMappings))

def _differentTranscript(prevAnnotRef, annot2GenomeRef):
    # ensure on same chrom for PAR issues
    return ((prevAnnotRef is not None) and
            ((annot2GenomeRef.transcriptId != prevAnnotRef.transcriptId) or
             (annot2GenomeRef.transcriptPos.name != prevAnnotRef.transcriptPos.name)))

def transAnnotMappingReader(annot2GenomePslFile, annot2GenomeRefTsv, annotLookupFunc,
                            transPslLookupFunc=None):
    """
    Reads mapped annotation alignments and metadata for target transcripts, including
    those that did not align successfully. Yields TransAnnotMapping objects.

    Args:
        annot2GenomePslFile: Path to a PSL file mapping annotations to the genome.
        annot2GenomeRefTsv: Path to a TSV file with annotation metadata.
        annotLookupFunc: Callable that takes an annotation ID and returns its data record.
        transPslLookupFunc: Optional callable that takes a transcript ID and chromosome,
            and returns the corresponding alignment information.

    Yields:
        TransAnnotMapping objects, one for each annotation.
    """

    annot2GenomePsls = [p for p in PslReader(annot2GenomePslFile)]
    prevAnnotRef = None
    annotMappings = []
    for annot2GenomeRef in annot2GenomeRefReader(annot2GenomeRefTsv):
        if _differentTranscript(prevAnnotRef, annot2GenomeRef):
            yield _makeTransAnnotMapping(annotMappings, transPslLookupFunc)
            annotMappings = []
        annotMappings.append(_makeAnnotMapping(annot2GenomeRef, annot2GenomePsls,
                                               annotLookupFunc(annot2GenomeRef.annotId)))
        prevAnnotRef = annot2GenomeRef
    if len(annotMappings) > 0:
        yield _makeTransAnnotMapping(annotMappings, transPslLookupFunc)
