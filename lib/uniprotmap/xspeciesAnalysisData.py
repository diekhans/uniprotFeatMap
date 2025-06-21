"""
Data structure to support xspecies analysis of annotations
"""
from dataclasses import dataclass
from collections import namedtuple, defaultdict
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
        # comment is allow to be None as wildcard
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


@dataclass
class SrcAnnotSet:
    """Uniprot source mapped annotations loaded into memory"""
    uniprotAnnotTbl: UniProtAnnotTbl
    annotMappingsTbl: AnnotMappingsTbl

def srcAnnotSetLoad(uniprotAnnotTsv, annot2GenomePsl, annot2GenomeRefTsv,
                    annotAssocs, targetGeneSet):
    def uniprotLookup(annotId):
        "None if annotation should be skipped"
        annot = uniprotAnnotTbl.getByAnnotId(annotId)
        return annot if annotAssocs.useSrc(annot.shortFeatType, annot.comment) else None

    uniprotAnnotTbl = UniProtAnnotTbl(uniprotAnnotTsv)
    annotMappingsTbl = transAnnotMappingLoader(annot2GenomePsl,
                                               annot2GenomeRefTsv,
                                               uniprotLookup,
                                               targetGeneSet.data.getAlign)
    return SrcAnnotSet(uniprotAnnotTbl, annotMappingsTbl)

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

    interproAnnotTbl = interproAnnotsLoad(interproAnnotTsv)
    annotMappingsTbl = transAnnotMappingLoader(interproAnnot2GenomePsl,
                                               interproAnnot2GenomeRefTsv,
                                               interproLookup,
                                               targetGeneSet.data.getAlign)
    return TargetAnnotSet(interproAnnotTbl, annotMappingsTbl)
