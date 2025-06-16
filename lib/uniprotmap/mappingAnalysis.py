"""
Analyze mappings of features to transcripts.
"""
from collections import namedtuple
from pycbio.sys.symEnum import SymEnum, auto

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

def _analyzeFullDeletion(annotRef):
    yield FeatureIndel(FeatureIndelType.del_full, annotRef.annotSize, 0, 0)

def _featureIndelGen(transPsl, annotMapping):
    if annotMapping.annotPsl is None:
        yield from _analyzeFullDeletion(annotMapping.annotRef)
    else:
        yield from _analyzePartialDeletion(transPsl, annotMapping.annotPsl)

def analyzeFeatureMapping(transAnnotMapping, annotMapping):
    """product list of feature disruptions, either
    unmapped regions of feature or insertions in
    feature that don't correspond to introns in
    transcripts."""
    return tuple(_featureIndelGen(transAnnotMapping.transPsl, annotMapping))
