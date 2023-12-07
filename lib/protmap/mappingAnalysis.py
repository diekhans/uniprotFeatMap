"""
Analyze mappings of features to transcripts.
"""
from collections import namedtuple
from pycbio.sys.symEnum import SymEnum, auto

class FeatureIndelType(SymEnum):
    del_5p = auto()
    del_3p = auto()
    del_int = auto()
    insert = auto()

_featIndelText = {
    FeatureIndelType.del_5p: "5' deletion",
    FeatureIndelType.del_3p: "5' deletion",
    FeatureIndelType.del_int: "internal deletion",
    FeatureIndelType.insert: "insertion"
}

def getFeatureIndelText(indelType):
    return _featIndelText[indelType]

class FeatureIndel(namedtuple("FeatureIndel",
                              ("indelType", "length",
                               "tStart", "tEnd",
                               ))):
    """Describes an INDEL.  For featDel, length is bases deleted, for inst it
    is length of insert.  ranges can be zero length."""
    __slots__ = ()

def _isIntron(transPsl, tStart, tEnd):
    "does the range match an intron in transPsl?"
    for iBlk in range(len(transPsl.blocks) - 1):
        iStart = transPsl.blocks[iBlk].tEnd
        iEnd = transPsl.blocks[iBlk + 1].tStart
        if (iStart == tStart) and (iEnd == tEnd):
            return True
        if iStart > tStart:
            break  # go past gap
    return False

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

def _analyzeBlock(transPsl, annotPsl, iBlk):
    blk = annotPsl.blocks[iBlk]
    nextBlk = annotPsl.blocks[iBlk + 1]
    if blk.qEnd < nextBlk.qStart:
        yield FeatureIndel(FeatureIndelType.del_int, (nextBlk.qStart - blk.qEnd),
                           blk.tEnd, nextBlk.tStart)
    if (blk.tEnd < nextBlk.tStart) and (not _isIntron(transPsl, blk.tEnd, nextBlk.tStart)):
        yield FeatureIndel(FeatureIndelType.insert, (nextBlk.tStart - blk.tEnd),
                           blk.tEnd, nextBlk.tStart)

def _analyzeBlocks(transPsl, annotPsl):
    for iBlk in range(len(annotPsl.blocks) - 1):
        yield from _analyzeBlock(transPsl, annotPsl, iBlk)

def _featureIndelGen(transPsl, annotPsl):
    yield from _analyzeStart(annotPsl)
    yield from _analyzeBlocks(transPsl, annotPsl)
    yield from _analyzeEnd(annotPsl)

def analyzeFeatureMapping(transPsl, annotPsl):
    """product list of feature disruptions, either
    unmapped regions of feature or insertions in
    feature that don't correspond to introns in
    transcripts."""
    return tuple(_featureIndelGen(transPsl, annotPsl))
