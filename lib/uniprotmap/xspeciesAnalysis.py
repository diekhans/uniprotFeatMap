"""
Compare mapped UniProt annotations to interpro annotations.
"""
from pycbio.sys.symEnum import SymEnum, auto
from uniprotmap.annotMappings import TransAnnotMappings
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
        """For deletions, coords covers range of deletion.
        """
        assert coords is not None
        # this the source, if there is one, otherwise there is a single target
        self.coords = coords
        # not initializes from arguments, as done in two phases
        self.srcAnnot = None
        self.srcFeatureIndels = None
        self.targetAnnots = []

    def dump(self, fh, transcriptId, indent=0):
        pre = indent * "  "
        pre2 = 2 * pre
        pre3 = 3 * pre
        print(f"{pre}diff: {transcriptId}: {self.coords}", file=fh)
        print(f"{pre2}srcAnnot: {self.srcAnnot}", file=fh)
        indelDesc = None
        if self.srcFeatureIndels is not None:
            indelDesc = [str(indel) for indel in self.srcFeatureIndels]
        print(f"{pre2}srcIndels: {indelDesc}", file=fh)
        print(f"{pre2}targetAnnot: {len(self.targetAnnots)}", file=fh)
        for targetAnnot in self.targetAnnots:
            print(f"{pre3}{targetAnnot}", file=fh)


class TransAnnotDiffs(list):
    """List of annotation comparisons differences for a transcript"""
    def __init__(self, transPsl):
        self.transPsl = transPsl

    def dump(self, fh, indent=0):
        print(f"TransAnnotDiffs: {self.transPsl.qName}", file=fh)
        for ad in self:
            ad.dump(fh, self.transPsl.qName, indent + 1)

def _makeSrcAnnotDiff(transAnnotMappings, annotMapping):
    """prevMappedTEnd & nextMappedTStart are use to fill in coords for region"""
    featureIndels = analyzeFeatureMapping(transAnnotMappings, annotMapping)
    diff = AnnotDiff(annotMapping.coords)
    diff.srcAnnot = annotMapping
    diff.srcFeatureIndels = featureIndels
    return diff

def _buildSrcDiffs(transAnnotDiffs, srcTransAnnotMappings):
    "must runs before adding target annots"
    for iAnnot, annotMapping in enumerate(srcTransAnnotMappings.annotMappings):
        transAnnotDiffs.append(_makeSrcAnnotDiff(srcTransAnnotMappings, annotMapping))

def _findTargetDiffIdx(transAnnotDiffs, annotMapping, iDiffs):
    """find point of first overlap in diffs lists, or before if no overlap """
    while iDiffs < len(transAnnotDiffs):
        if annotMapping.coords.start < transAnnotDiffs[iDiffs].coords.end:
            if transAnnotDiffs[iDiffs].srcAnnot is not None:
                return iDiffs
        iDiffs += 1
    return iDiffs  # add to end

def _targetDiffNew(annotMapping):
    diff = AnnotDiff(annotMapping.coords)
    diff.targetAnnots.append(annotMapping)
    return diff

def _addTargetDiffOverlap(transAnnotDiffs, annotMapping, iDiffs):
    while iDiffs < len(transAnnotDiffs) and annotMapping.coords.overlaps(transAnnotDiffs[iDiffs].coords):
        assert transAnnotDiffs[iDiffs].srcAnnot is not None
        transAnnotDiffs[iDiffs].targetAnnots.append(annotMapping)
        iDiffs += 1

def _addTargetDiff(transAnnotDiffs, annotMapping, iDiffs):
    iDiffs = _findTargetDiffIdx(transAnnotDiffs, annotMapping, iDiffs)
    if iDiffs >= len(transAnnotDiffs):
        transAnnotDiffs.targetAnnots.append(_targetDiffNew(annotMapping))
    elif not annotMapping.coords.overlaps(transAnnotDiffs[iDiffs].coords):
        transAnnotDiffs.insert(iDiffs, _targetDiffNew(annotMapping))
    else:
        _addTargetDiffOverlap(transAnnotDiffs, annotMapping, iDiffs)
    return iDiffs

def _addTargetDiffs(transAnnotDiffs, targetTransAnnotMappings):
    "add in the target (interpro) annotations"
    iDiffs = 0  # point to start search for overlap or insert
    for annotMapping in targetTransAnnotMappings.annotMappings:
        iDiffs = _addTargetDiff(transAnnotDiffs, annotMapping, iDiffs)

def compareTransAnnotations(srcTransAnnotMappings: TransAnnotMappings,
                            targetTransAnnotMappings: TransAnnotMappings) -> TransAnnotDiffs:
    """build comparison of src to target domain annotations.  Inputs must be sorted"""
    transAnnotDiffs = TransAnnotDiffs(srcTransAnnotMappings.transPsl)
    _buildSrcDiffs(transAnnotDiffs, srcTransAnnotMappings)
    _addTargetDiffs(transAnnotDiffs, targetTransAnnotMappings)
    return transAnnotDiffs
