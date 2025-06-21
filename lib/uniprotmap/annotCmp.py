"""
Compares annotations on transcripts from different sources.
"""
from dataclasses import dataclass
from pycbio.hgdata.coords import Coords
from uniprotmap.geneset import geneSetFactory
from uniprotmap.uniprot import UniProtAnnotTbl
from uniprotmap.interproscan import InterproAnnotTbl, interproAnnotsLoad
from uniprotmap.annotMappings import AnnotMappingsTbl, transAnnotMappingLoader
from uniprotmap.mappingAnalysis import analyzeFeatureMapping

class AnnotAssoc(namedtuple("AnnotAssoc",
                             ("uniprotShortFeatType", "interproAnalysis", "interproAcc"))):
    """
    Specifies a Uniprot anotation and the match interpro analysis/acc.
    """

class AnnotAssocs:
    """Takes list of AnnotAssoc objects which are sets that are considered
    equivalent.  UniProt annotations maybe not be overlapping. Multiple
    interpro annotations with the same accession from different analysis may
    overlap.
    """

    def __init__(self):
        self.srcToTargets = defaultdict(list)
        self.targetToSrcs = defaultdict(list)  # by (analysis, acc)

    def add(self, uniprotShortFeatType, interproAnalysis, interproAcc):
        annotAssoc = AnnotAssoc(uniprotShortFeatType, interproAnalysis, interproAcc)
        self.srcToTargets[uniprotShortFeatType].add(annotAssoc)
        self.targetToSrcs[(interproAnalysis, interproAcc)].add(annotAssoc)

    def finish(self):
        self.srcToTargets.default_factory = self.targetToSrcs.default_factory = None

    def useSrc(self, uniprotShortFeatType):
        return uniprotShortFeatType in self.srcToTargets

    def useTarget(self, interproAnalysis, interproAcc):
        return (interproAnalysis, interproAcc) in self.targetToSrcs

    def getSrcTargets(self, uniprotShortFeatType):
        return self.srcToTargets[uniprotShortFeatType]

@dataclass
class SrcAnnotSet:
    """Uniprot source mapped annotations loaded into memory"""
    uniprotAnnotTbl: UniProtAnnotTbl
    annotMappingsTbl: AnnotMappingsTbl

def srcAnnotSetLoad(srcAnnotConf, annotFilters, targetGeneSet):
    def uniprotLookup(annotId):
        annot = uniprotAnnotTbl.getByAnnotId(annotId)
        return annot if annotFilters.useSrc(annot.shortFeatType) else None

    def pslLookup(transId, chrom):
        return targetGeneSet.getEntry(transId, chrom).psl

    uniprotAnnotTbl = UniProtAnnotTbl(srcAnnotConf.uniprotAnnotTsv)
    annotMappingsTbl = transAnnotMappingLoader(srcAnnotConf.annot2GenomePsl,
                                               srcAnnotConf.annot2GenomeRefTsv,
                                               uniprotLookup, pslLookup)
    return SrcAnnotSet(uniprotAnnotTbl, annotMappingsTbl)

@dataclass
class TargetAnnotSet:
    """Set to evaluate against source"""
    annotData: InterproAnnotTbl
    annotMappingsTbl: AnnotMappingsTbl

def targetAnnotSetLoad(targetAnnotConf):
    annotData = interproAnnotsLoad(targetAnnotConf.annotTsv)
    return TargetAnnotSet(annotData)


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
    for annotMapping in transAnnotMappings.annotMappings:
        pass
