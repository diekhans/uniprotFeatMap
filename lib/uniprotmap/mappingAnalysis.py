"""
Analyze mappings of features to transcripts.
"""
from collections import namedtuple, defaultdict
from pycbio.sys.symEnum import SymEnum, auto
from pycbio.hgdata.psl import PslReader
from uniprotmap.metadata import annot2GenomeRefReader

class MappingError(Exception):
    pass


###
# Reading annotation mappings
###
class AnnotMapping(namedtuple("AnnotMapping",
                              ("annotRef", "annotPsl", "annot"))):
    """A single mapping for an annotation to a genome. annotPsl is None if not mapped.
    The annot field depends on the type of annotation (UniProt or Interpro).
    """
    __slots__ = ()


class TransAnnotMappings(namedtuple("TransAnnotMappings",
                                    ("transcriptId", "chrom", "annotMappings",))):
    """Mappings for all annotations on a transcript.  The full list is
    needed to place deleted annotations."""
    __slots__ = ()

class AnnotMappingsTbl(list):
    """Table of TransAnnotMappings, all the mappings annotations to a genome"""
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


def _makeAnnotMapping(annot2GenomeRef, annot2GenomePsls, annot):
    if annot2GenomeRef.alignIdx is None:
        return AnnotMapping(annot2GenomeRef, None, annot)
    else:
        return AnnotMapping(annot2GenomeRef, annot2GenomePsls[annot2GenomeRef.alignIdx], annot)

def _makeTransAnnotMapping(annotMappings):
    annotRef = annotMappings[0].annotRef
    return TransAnnotMappings(annotRef.transcriptId,
                              annotRef.transcriptPos.name,
                              tuple(annotMappings))

def _differentTranscript(prevAnnotRef, annot2GenomeRef):
    # ensure on same chrom for PAR issues
    return ((prevAnnotRef is not None) and
            ((annot2GenomeRef.transcriptId != prevAnnotRef.transcriptId) or
             (annot2GenomeRef.transcriptPos.name != prevAnnotRef.transcriptPos.name)))
def transAnnotMappingReader(annot2GenomePslFile, annot2GenomeRefTsv, annotTbl):
    """Reads mapped annotation alignments and metadata for a target transcript.  Returns only
    metadata for annotations that don't map. Yields TransAnnotMappings objects,
    The annotTbl.getByAnnotId() function is given the annotation id and should return the specific
    annotation data.
    """

    annot2GenomePsls = [p for p in PslReader(annot2GenomePslFile)]
    prevAnnotRef = None
    annotMappings = []
    for annot2GenomeRef in annot2GenomeRefReader(annot2GenomeRefTsv):
        if _differentTranscript(prevAnnotRef, annot2GenomeRef):
            yield _makeTransAnnotMapping(annotMappings)
            annotMappings = []
        annotMappings.append(_makeAnnotMapping(annot2GenomeRef, annot2GenomePsls,
                                               annotTbl.getByAnnotId(annot2GenomeRef.annotId)))
        prevAnnotRef = annot2GenomeRef
    if len(annotMappings) > 0:
        yield _makeTransAnnotMapping(annotMappings)


###
# Analysis of annotation INDELs
###
class FeatureIndelType(SymEnum):
    del_5p = auto()
    del_3p = auto()
    del_int = auto()
    del_full = auto()
    insert = auto()

_featIndelText = {
    FeatureIndelType.del_5p: "5' deletion",
    FeatureIndelType.del_3p: "5' deletion",
    FeatureIndelType.del_int: "internal deletion",
    FeatureIndelType.del_full: "full deletion",
    FeatureIndelType.insert: "insertion"
}

def getFeatureIndelText(indelType):
    return _featIndelText[indelType]

class FeatureIndel(namedtuple("FeatureIndel",
                              ("indelType", "length", "tStart", "tEnd"))):
    """Describes an INDEL.  For featDel, length in bases deleted, for insert it
    is length of insert.  ranges can be zero length."""
    __slots__ = ()

def _exonIntersections(transPsl, tStart, tEnd):
    """return list of genome ranges where exons intersect the specified range.
    If nothing is return, that indicates tStart/tEnd are exactly an exon.
    """
    exonIntersects = []
    for iBlk in range(len(transPsl.blocks)):
        start = max(tStart, transPsl.blocks[iBlk].tStart)
        end = min(tEnd, transPsl.blocks[iBlk].tEnd)
        if start < end:
            exonIntersects.append((start, end))
    return exonIntersects

def _analyzeStart(annotPsl):
    if annotPsl.qStart > 0:
        # unaligned before start of query
        if annotPsl.qStrand == '+':
            tPos = annotPsl.blocks[0].tStart
            indelType = FeatureIndelType.del_5p
        else:
            tPos = annotPsl.blocks[-1].tEnd
            indelType = FeatureIndelType.del_3p
        yield FeatureIndel(indelType, annotPsl.qStart, tPos, tPos)

def _analyzeEnd(annotPsl):
    if annotPsl.qEnd < annotPsl.qSize:
        # unaligned after last block
        if annotPsl.qStrand == '+':
            tPos = annotPsl.blocks[-1].tEnd
            indelType = FeatureIndelType.del_3p
        else:
            tPos = annotPsl.blocks[0].tStart
            indelType = FeatureIndelType.del_5p
        yield FeatureIndel(indelType, annotPsl.qSize - annotPsl.qEnd, tPos, tPos)

def _analyzeQDel(blk, nextBlk):
    yield FeatureIndel(FeatureIndelType.del_int, (nextBlk.qStart - blk.qEnd),
                       blk.tEnd, nextBlk.tStart)

def _analyzeTDel(transPsl, blk, nextBlk):
    # this excludes introns
    for exStart, exEnd in _exonIntersections(transPsl, blk.tEnd, nextBlk.tStart):
        yield FeatureIndel(FeatureIndelType.insert, (exEnd - exStart),
                           exStart, exEnd)

def _analyzeBlock(transPsl, annotPsl, iBlk):
    blk = annotPsl.blocks[iBlk]
    nextBlk = annotPsl.blocks[iBlk + 1]
    if blk.qEnd < nextBlk.qStart:
        yield from _analyzeQDel(blk, nextBlk)

    if blk.tEnd < nextBlk.tStart:
        yield from _analyzeTDel(transPsl, blk, nextBlk)

def _analyzeBlocks(transPsl, annotPsl):
    for iBlk in range(len(annotPsl.blocks) - 1):
        yield from _analyzeBlock(transPsl, annotPsl, iBlk)

def _analyzePartialDeletion(transPsl, annotPsl):
    yield from _analyzeStart(annotPsl)
    yield from _analyzeBlocks(transPsl, annotPsl)
    yield from _analyzeEnd(annotPsl)

def _analyzeFullDeletion():
    raise Exception("_analyzeFullDeletion not implemented yet")

def _featureIndelGen(transPsl, annotMapping):
    if annotMapping.annotPsl is None:
        yield from _analyzeFullDeletion()
    else:
        yield from _analyzePartialDeletion(transPsl, annotMapping.annotPsl)

def analyzeFeatureMapping(transPsl, annotMapping):
    """product list of feature disruptions, either
    unmapped regions of feature or insertions in
    feature that don't correspond to introns in
    transcripts."""
    return tuple(_featureIndelGen(transPsl, annotMapping))
