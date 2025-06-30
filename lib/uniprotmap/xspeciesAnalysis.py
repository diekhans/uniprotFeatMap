"""
Compare mapped UniProt annotations to interpro annotations.
"""
from pycbio.hgdata.coords import Coords
from pycbio.sys.symEnum import SymEnum, auto
from uniprotmap.mappingAnalysis import analyzeFeatureMapping

class AnnotMethod(SymEnum):
    """how was annotation being compared produced?"""
    uniprotMap = auto()
    interpro = auto()

class AnnotDiffCategory(SymEnum):
    """Summary of the AnnotDiff, not all valid for both src and target"""
    complete = auto()
    deleted = auto()
    inserted = auto()   # for inserts by target
    minor_diff = auto()
    major_diff = auto()


class AnnotDiff:
    """describes diffs this can be one src and overlapping targets """
    def __init__(self, coords):
        """For deletions, coords covers range of deletion
        """
        self.coords = coords
        # not initializes from arguments, as done in two phases
        self.srcAnnot = None
        self.srcFeatureIndels = None
        self.targetAnnot = None

    def dump(self, fh, transcriptId, indent=0):
        pre = indent * "  "
        pre2 = 2 * pre
        print(f"{pre}diff: {transcriptId}: {self.coords}", file=fh)
        print(f"{pre2}srcAnnot: {self.srcAnnot}", file=fh)
        indelDesc = None
        if self.srcFeatureIndels is not None:
            indelDesc = [str(indel) for indel in self.srcFeatureIndels]
        print(f"{pre2}srcIndels: {indelDesc}", file=fh)
        print(f"{pre2}targetAnnot: {self.targetAnnot}", file=fh)

class TransAnnotDiffs(list):
    """Annotation comparisons differences for a transcript"""
    def __init__(self, transPsl):
        self.transPsl = transPsl

    def dump(self, fh, indent=0):
        print(f"TransAnnotDiffs: {self.transPsl.qName}", file=fh)
        for ad in self:
            ad.dump(fh, self.transPsl.qName, indent + 1)

def _makeSrcAnnotDiff(transAnnotMappings, annotMapping, prevMappedTEnd, nextMappedTStart):
    """prevMappedTEnd & nextMappedTStart are use to fill in coords for region"""
    featureIndels = analyzeFeatureMapping(transAnnotMappings, annotMapping)
    if annotMapping.annotPsl is not None:
        coords = Coords(annotMapping.annotPsl.tName, annotMapping.annotPsl.tStart, annotMapping.annotPsl.tEnd)
    else:
        # range that could contaion features
        coords = Coords(transAnnotMappings.transPsl.tName, prevMappedTEnd, nextMappedTStart)
    entry = AnnotDiff(coords)
    entry.srcAnnot = annotMapping
    entry.srcFeatureIndels = featureIndels
    return entry

def _findNextMapped(transAnnotMappings, iAnnot):
    """next non-deletion, returning start coords, or end of transcript in no more"""
    iAnnot += 1
    while iAnnot < len(transAnnotMappings.annotMappings):
        annotMapping = transAnnotMappings.annotMappings[iAnnot]
        if annotMapping.annotPsl is not None:
            return annotMapping.annotPsl.tStart
        iAnnot += 1
    return transAnnotMappings.transPsl.tEnd

def _buildSrcDiffs(transAnnotDiffs, srcTransAnnotMappings):
    "runs before adding target annots"
    # build initial list of entries, prev/next is use to find location
    # where a deleted could be
    prevMappedTEnd = srcTransAnnotMappings.transPsl.tStart
    for iAnnot, annotMapping in enumerate(srcTransAnnotMappings.annotMappings):
        entry = _makeSrcAnnotDiff(srcTransAnnotMappings, annotMapping, prevMappedTEnd,
                                  _findNextMapped(srcTransAnnotMappings, iAnnot))
        transAnnotDiffs.append(entry)
        if annotMapping.annotPsl is not None:
            prevMappedTEnd = annotMapping.annotPsl.tEnd

def _findTargetDiffIdx(transAnnotDiffs, annotMapping, iDiffs):
    """find point of first overlap in diffs lists, or before if no overlap """
    annotCoords = annotMapping.coords
    while iDiffs < len(transAnnotDiffs):
        pass


def _addTargetDiff(transAnnotDiffs, annotMapping, iDiffs):
    pass


def _addTargetDiffs(transAnnotDiffs, targetTransAnnotMappings):
    "add in the target (interpro) annotations"
    iDiffs = 0  # point to start search for overlap or insert
    for annotMapping in targetTransAnnotMappings.annotMappings:
        pass
        iDiffs = _addTargetDiff(transAnnotDiffs, annotMapping, iDiffs)


def compareTransAnnotations(srcTransAnnotMappings,
                            targetTransAnnotMappings):
    transAnnotDiffs = TransAnnotDiffs(srcTransAnnotMappings.transPsl)
    _buildSrcDiffs(transAnnotDiffs, srcTransAnnotMappings)
    return transAnnotDiffs
