"""
Analyze mappings of features to transcripts.
"""
from collections import namedtuple

class FeatureIndel(namedtuple("FeatureIndel",
                              ("isFeatDel", "length",
                               "qStart", "qEnd",
                               "tStart", "tEnd",
                               ))):
    """Describes an INDEL.  For featDel, length is
    bases deleted, for inst it is length of insert.
    ranges can be zero length"""
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
    if annotPsl.blocks[0].qStart > 0:
        yield FeatureIndel(True, annotPsl.blocks[0].qStart,#FIXME: wrong, needs to be realtive to qStart
                           0, annotPsl.blocks[0].qStart,
                           annotPsl.blocks[0].tStart, annotPsl.blocks[0].tStart)

def _analyzeEnd(annotPsl):
    if annotPsl.blocks[-1].qEnd < annotPsl.qSize:
        yield FeatureIndel(True, annotPsl.qSize - annotPsl.blocks[-1].qEnd,#WRONG if strand
                           annotPsl.blocks[-1].qEnd, annotPsl.qSize,
                           annotPsl.blocks[-1].tEnd, annotPsl.blocks[-1].tEnd)

def _analyzeBlock(transPsl, annotPsl, iBlk):
    blk = annotPsl.blocks[iBlk]
    nextBlk = annotPsl.blocks[iBlk + 1]
    if blk.qEnd < nextBlk.qStart:
        yield FeatureIndel(True, (nextBlk.qStart - blk.qEnd),
                           blk.qEnd, nextBlk.qStart,
                           blk.tEnd, nextBlk.tStart)
    if (blk.tEnd < nextBlk.tStart) and (not _isIntron(transPsl, blk.tEnd, nextBlk.tStart)):
        yield FeatureIndel(False, (nextBlk.tStart - blk.tEnd),
                           blk.qEnd, nextBlk.qStart,
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
