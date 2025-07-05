"""
Data structure to support xspecies analysis of annotations
"""
from dataclasses import dataclass
from collections import namedtuple, defaultdict
from itertools import islice
from uniprotmap import DataError
from uniprotmap.uniprot import UniProtAnnotTbl
from uniprotmap.interproscan import InterproAnnotTbl, interproAnnotsLoad
from uniprotmap.annotMappings import AnnotMappingsTbl, transAnnotMappingLoader

class AnnotAssoc(namedtuple("AnnotAssoc",
                            ("uniprotShortFeatType", "uniprotComment",
                             "interproAnalysis", "interproAcc"))):
    """
    Specifies a Uniprot anotation and the match interpro analysis/acc.
    Uniprot comment is needed to distinguish the type of some domains.
    If uniprotComment is None, it is not matched.
    """
    pass

class AnnotAssocs:
    """Takes list of AnnotAssoc objects which are sets that are considered
    equivalent.  UniProt annotations maybe not be overlapping. Multiple
    interpro annotations with the same accession from different analysis may
    overlap.
    """

    def __init__(self):
        self.srcToTargets = defaultdict(list)  # by (shortType, comment
        self.targetToSrcs = defaultdict(list)  # by (analysis, acc)

    def add(self, uniprotShortFeatType, uniprotComment, interproAnalysis, interproAcc):
        annotAssoc = AnnotAssoc(uniprotShortFeatType, uniprotComment,
                                interproAnalysis, interproAcc)
        self.srcToTargets[(uniprotShortFeatType, uniprotComment)].append(annotAssoc)
        self.targetToSrcs[(interproAnalysis, interproAcc)].append(annotAssoc)

    def finish(self):
        self.srcToTargets.default_factory = self.targetToSrcs.default_factory = None

    def useSrc(self, uniprotShortFeatType, uniprotComment):
        # None comment is a while card
        return (((uniprotShortFeatType, uniprotComment) in self.srcToTargets) or
                ((uniprotShortFeatType, None) in self.srcToTargets))

    def useTarget(self, interproAnalysis, interproAcc):
        return (interproAnalysis, interproAcc) in self.targetToSrcs

    def getSrcTargets(self, uniprotShortFeatType, uniprotComment):
        targets = self.srcToTargets.get((uniprotShortFeatType, uniprotComment))
        if targets is None:
            targets = self.srcToTargets.get((uniprotShortFeatType, None))
        assert targets is not None
        return targets

def _raiseOverlapError(prevMapping, annotMapping, msg):
    raise DataError(f"{msg}:\n`" + str(prevMapping), "', and `" + str(annotMapping) + "'")

def _checkForOverlapTransAnnot(prevMapping, annotMapping):
    if annotMapping.coords.overlaps(prevMapping.coords):
        _raiseOverlapError(prevMapping, annotMapping,
                           "overlapping annotations not allowed, should be filtered")
    if annotMapping.coords.start < prevMapping.coords.end:
        _raiseOverlapError(prevMapping, annotMapping,
                           "out of order annotations")

def _checkForOverlapTransAnnots(transAnnotMapping):
    prevMapping = transAnnotMapping.annotMappings[0]
    for annotMapping in islice(transAnnotMapping.annotMappings, 1, None):
        if (prevMapping.coords is not None) and (annotMapping.coords is not None):
            _checkForOverlapTransAnnot(prevMapping, annotMapping)
        prevMapping = annotMapping

def _checkForOverlapAnnots(annotMappingsTbl):
    """code doesn't handle overlapping source mapped annotations on a transcript at
    this point"""
    for transAnnotMapping in annotMappingsTbl:
        _checkForOverlapTransAnnots(transAnnotMapping)

@dataclass
class SrcAnnotSet:
    """Uniprot source mapped annotations loaded into memory"""
    uniprotAnnotTbl: UniProtAnnotTbl
    annotMappingsTbl: AnnotMappingsTbl

def srcAnnotSetLoad(uniprotAnnotTsv, uniprotAnnot2GenomePsl, uniprotAnnot2GenomeRefTsv,
                    annotAssocs, targetGeneSet):
    def uniprotLookup(annotId):
        "None if annotation should be skipped"
        annot = uniprotAnnotTbl.getByAnnotId(annotId)
        return annot if annotAssocs.useSrc(annot.shortFeatType, annot.comment) else None

    try:
        uniprotAnnotTbl = UniProtAnnotTbl(uniprotAnnotTsv)
        annotMappingsTbl = transAnnotMappingLoader(uniprotAnnot2GenomePsl,
                                                   uniprotAnnot2GenomeRefTsv,
                                                   uniprotLookup,
                                                   targetGeneSet.data.getAlign)
        _checkForOverlapAnnots(annotMappingsTbl)
        return SrcAnnotSet(uniprotAnnotTbl, annotMappingsTbl)
    except Exception as ex:
        raise DataError("problem load mapped source annotations from `" +
                        uniprotAnnotTsv + "', `" +
                        uniprotAnnot2GenomePsl + "', and `" +
                        uniprotAnnot2GenomeRefTsv + "'") from ex

@dataclass
class TargetAnnotSet:
    """Set to evaluate against source"""
    intreproAnnotTbl: InterproAnnotTbl
    annotMappingsTbl: AnnotMappingsTbl

def targetAnnotSetLoad(interproAnnotTsv, interproAnnot2GenomePsl, interproAnnot2GenomeRefTsv,
                       annotAssocs, targetGeneSet):
    def interproLookup(annotId):
        "None if annotation should be skipped"
        annot = interproAnnotTbl.getByAnnotId(annotId)
        return annot if annotAssocs.useTarget(annot.analysis, annot.interpro_accession) else None
    try:
        interproAnnotTbl = interproAnnotsLoad(interproAnnotTsv)
        annotMappingsTbl = transAnnotMappingLoader(interproAnnot2GenomePsl,
                                                   interproAnnot2GenomeRefTsv,
                                                   interproLookup, targetGeneSet.data.getAlign,
                                                   sortByCoords=True)
        return TargetAnnotSet(interproAnnotTbl, annotMappingsTbl)
    except Exception as ex:
        raise DataError("problem load target annotations from `" +
                        interproAnnotTsv + "', `" +
                        interproAnnot2GenomePsl + "', and `" +
                        interproAnnot2GenomeRefTsv + "'") from ex
