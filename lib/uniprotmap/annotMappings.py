"""
Annotation mappings data for analysis
"""
from collections import namedtuple, defaultdict
from pycbio.hgdata.psl import PslReader
from pycbio.hgdata.coords import Coords
from uniprotmap.metadata import annot2GenomeRefReader


class MappingError(Exception):
    pass


class AnnotMapping(namedtuple("AnnotMapping",
                              ("annotRef", "annotPsl", "annot", "coords"))):
    """A single mapping for an annotation to a genome. annotPsl is None if not
    mapped.  The annot field depends on the type of annotation (UniProt or
    Interpro).  The coords field is either the range from the PSL, or if not
    mapped, the gap between mapped annotations or the ends of the transcript,
    where the annotations would have mapped.
    """
    __slots__ = ()

    def __new__(cls, annotRef, annotPsl, annot, coords=None):
        """Option coords will be filled in from PSL if not specified.  It is needed
        as an argument for unpickle"""
        if (coords is None) and (annotPsl is not None):
            coords = Coords(annotPsl.tName, annotPsl.tStart, annotPsl.tEnd,
                            annotPsl.tStrand, annotPsl.tSize)
        return super(AnnotMapping, cls).__new__(cls, annotRef, annotPsl, annot, coords)

    def short(self):
        return f"{self.annotRef.annotId}/{self.annotRef.annotId}/{self.annotRef.annotMapId}: {self.annot.short()} {self.coords}"


class TransAnnotMappings(namedtuple("TransAnnotMappings",
                                    ("transcriptId", "chrom", "transPsl", "annotMappings",))):
    """Mappings for all annotations on a transcript."""
    __slots__ = ()


class AnnotMappingsTbl(list):
    """Table of TransAnnotMappings, all the mappings annotations to a genome"""
    def __init__(self):
        self.byTransId = defaultdict(list)   # multiple to handle PAR

    def add(self, transAnnotMappings):
        self.append(transAnnotMappings)
        self.byTransId[transAnnotMappings.transcriptId].append(transAnnotMappings)

    def finish(self):
        self.byTransId.default_factory = None

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
def _makeAnnotMapping(annot2GenomeRef, annot2GenomePsls, annotLookupFunc):
    annotPsl = coords = None
    if annot2GenomeRef.alignIdx is not None:
        annotPsl = annot2GenomePsls[annot2GenomeRef.alignIdx]
        coords = Coords(annotPsl.tName, annotPsl.tStart, annotPsl.tEnd)

    # this might filter annotation
    annot = annotLookupFunc(annot2GenomeRef.annotId)
    if annot is not None:
        return AnnotMapping(annot2GenomeRef, annotPsl, annot, coords)
    else:
        return None

def _makeAnnotMappings(transAnnot2GenomeRefs, annot2GenomePsls, annotLookupFunc):
    annotMappings = []
    for annot2GenomeRef in transAnnot2GenomeRefs:
        annotMapping = _makeAnnotMapping(annot2GenomeRef, annot2GenomePsls, annotLookupFunc)
        if annotMapping is not None:
            annotMappings.append(annotMapping)
    return annotMappings

def _makeTransAnnotMapping(transAnnot2GenomeRefs, annot2GenomePsls, annotLookupFunc, transPslLookupFunc,
                           inTranscriptionOrder):
    annotRef0 = transAnnot2GenomeRefs[0]
    transPsl = transPslLookupFunc(annotRef0.transcriptId,
                                  annotRef0.transcriptPos.name)
    annotMappings = _makeAnnotMappings(transAnnot2GenomeRefs, annot2GenomePsls, annotLookupFunc)
    if inTranscriptionOrder and (transPsl.qStrand == '-'):
        annotMappings.reverse()
    return TransAnnotMappings(annotRef0.transcriptId,
                              annotRef0.transcriptPos.name,
                              transPsl, tuple(annotMappings))

def _differentTranscript(prevAnnotRef, annot2GenomeRef):
    # ensure on same chrom for PAR issues
    return ((prevAnnotRef is not None) and
            ((annot2GenomeRef.transcriptId != prevAnnotRef.transcriptId) or
             (annot2GenomeRef.transcriptPos.name != prevAnnotRef.transcriptPos.name)))

def _transAnnot2GenomeRefReader(annot2GenomeRefTsv):
    """returns list of annotations associated with the next transcript"""
    prevAnnotRef = None
    transAnnot2GenomeRefs = []
    for annot2GenomeRef in annot2GenomeRefReader(annot2GenomeRefTsv):
        if _differentTranscript(prevAnnotRef, annot2GenomeRef):
            yield transAnnot2GenomeRefs
            transAnnot2GenomeRefs = []
        transAnnot2GenomeRefs.append(annot2GenomeRef)
        prevAnnotRef = annot2GenomeRef
    if len(transAnnot2GenomeRefs) > 0:
        yield transAnnot2GenomeRefs

def transAnnotMappingReader(annot2GenomePslFile, annot2GenomeRefTsv, annotLookupFunc,
                            transPslLookupFunc, *, inTranscriptionOrder=False):
    """
    Reads mapped annotation alignments and metadata for target transcripts, including
    those that did not align successfully. Yields TransAnnotMapping objects.

    Args:
        annot2GenomePslFile: Path to a PSL file mapping annotations to the genome.
        annot2GenomeRefTsv: Path to a TSV file with annotation metadata.
        annotLookupFunc: Callable that takes an annotation ID and returns its data record.
            It is also a filter, if it returns None, the annotation mapping is discarded.
        transPslLookupFunc: Callable that takes a transcript ID and chromosome,
            and returns the corresponding alignment information.
        sortByCoords: If True. then sort by coordinates.  Can not be used on spared mapped
            annotations.
    Yields:
        TransAnnotMapping objects, one for each annotation.
    """

    annot2GenomePsls = [p for p in PslReader(annot2GenomePslFile)]
    for transAnnot2GenomeRefs in _transAnnot2GenomeRefReader(annot2GenomeRefTsv):
        yield _makeTransAnnotMapping(transAnnot2GenomeRefs, annot2GenomePsls, annotLookupFunc, transPslLookupFunc,
                                     inTranscriptionOrder)


def transAnnotMappingLoader(annot2GenomePslFile, annot2GenomeRefTsv, annotLookupFunc,
                            transPslLookupFunc, *, inTranscriptionOrder=False):
    """load mappings into object AnnotMappingsTbl"""
    annotMappingsTbl = AnnotMappingsTbl()
    for transAnnotMappings in transAnnotMappingReader(annot2GenomePslFile, annot2GenomeRefTsv,
                                                      annotLookupFunc, transPslLookupFunc,
                                                      inTranscriptionOrder=inTranscriptionOrder):
        annotMappingsTbl.add(transAnnotMappings)
    annotMappingsTbl.finish()
    return annotMappingsTbl
