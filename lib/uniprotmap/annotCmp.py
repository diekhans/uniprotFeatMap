"""
Compares annotations on transcripts from different sources.
"""
from pycbio.hgdata.coords import Coords
from uniprotmap.mappingAnalysis import analyzeFeatureMapping

class AnnotCmpEntry:
    """describes set of mappings"""
    def __init__(self, coords):
        "for deletions, covers possible locations for deletion"
        self.coords = coords
        self.srcAnnot = None
        self.srcFeatureIndels = None
        self.targetAnnots = None

class TransAnnotCmp:
    """Annotation comparisons for a transcript"""
    def __init__(self, transPsl, entries):
        self.transPsl = transPsl
        self.entries = entries


def _makeSrcAnnotCmp(transAnnotMappings, annotMapping, prevMappedTEnd, nextMappedTStart):
    """prevMappedTEnd & nextMappedTStart are use to fill in coords for region"""
    featureIndels = None
    if annotMapping.annotPsl is not None:
        coords = Coords(annotMapping.annotPsl.tName, annotMapping.annotPsl.tStart, annotMapping.annotPsl.tEnd)
        featureIndels = analyzeFeatureMapping(transAnnotMappings, annotMapping)
    else:
        coords = Coords(transAnnotMappings.transPsl.tName, prevMappedTEnd, nextMappedTStart)
    entry = AnnotCmpEntry(coords)
    entry.srcAnnot = annotMapping
    entry.srcFeatureIndels = featureIndels
    return entry

def _findNextMapped(transAnnotMappings, iAnnot):
    """next non-deletion, returning start coords, or end of transcript in no more"""
    iAnnot += 1
    while iAnnot < transAnnotMappings.annotMappings:
        annotMapping = transAnnotMappings.annotMappings[iAnnot]
        if annotMapping.annotPsl is not None:
            return annotMapping.annotPsl.tStart
    return transAnnotMappings.transPsl.tEnd

def _buildSrcAnnotMappingsCmp(transAnnotMappings):
    "runs before adding target annots"
    # build initial list of entries, prev/next is use to find location
    # where a deleted could be
    prevMappedTEnd = transAnnotMappings.transPsl.tStart
    entries = []
    for iAnnot, annotMapping in enumerate(transAnnotMappings.annotMappings):
        entry = _makeSrcAnnotCmp(transAnnotMappings, annotMapping, prevMappedTEnd,
                                 _findNextMapped(transAnnotMappings, iAnnot))
        entries.append(entry)
        if annotMapping.annotPsl is not None:
            prevMappedTEnd = annotMapping.annotPsl.tEnd
    return TransAnnotCmp(transAnnotMappings.transPsl, entries)

def _addTargetAnnotsCmp(targetTransAnnotMappings, transAnnotCmp):
    "add in the target (interpro) annotations"
    iNext = 0
    for annotMapping in targetTransAnnotMappings.annotMappings:
        pass
