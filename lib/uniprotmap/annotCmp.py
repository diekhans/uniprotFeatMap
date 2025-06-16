"""
Compares annotations on transcripts from different sources.
"""
from dataclasses import dataclass
from typing import Union
from pycbio.hgdata.psl import PslTbl
from uniprotmap.annotCmpConf import AnnotSet
from uniprotmap.uniprot import UniProtAnnotTbl
from uniprotmap.interproscan import InterproAnnotTbl, interproAnnotsLoad
from uniprotmap.mappingAnalysis import AnnotMappingsTbl

@dataclass
class SrcAnnotSet:
    """Source annotations loaded into memory"""
    annotSet: AnnotSet
    uniprotAnnotTbl: UniProtAnnotTbl
    annotMappingsTbl: AnnotMappingsTbl
    xspeciesTrans2TransTbl: PslTbl

def srcAnnotSetLoad(srcAnnotConf):
    return SrcAnnotSet(srcAnnotConf.annotSet,
                       UniProtAnnotTbl(srcAnnotConf.uniprotAnnotTsv),
                       AnnotMappingsTbl(srcAnnotConf.annot2GenomePsl,
                                        srcAnnotConf.annot2GenomeRefTsv),
                       PslTbl(srcAnnotConf.xspeciesTrans2TransPsl, qNameIdx=True))

@dataclass
class TargetAnnotSet:
    """Set to evaluate against source"""
    annotSet: AnnotSet
    annotData: Union[UniProtAnnotTbl, InterproAnnotTbl]

def targetAnnotSetLoad(targetAnnotConf):
    if targetAnnotConf.annotSet is AnnotSet.InterPro:
        annotData = UniProtAnnotTbl(targetAnnotConf.annotTsv)
    else:
        annotData = interproAnnotsLoad(targetAnnotConf.annotTsv)
    return TargetAnnotSet(targetAnnotConf.annotSet, annotData)

# def

    # annotTsv: str    # either uniprotAnnotTsv or interproAnnotTsv output
    # annot2TransPsl: str
    # annot2GenomePsl: str
    # annot2GenomeRefTsv: str


class AnnotMap:
    """Mappings of an annotation to target genome"""
    __slots__ = ('annotId', 'mapRanges')

    def __init__(self, annotId, mapRanges):
        self.annotId = annotId
        self.mapRanges = mapRanges

class SrcTranscript:
    "source transcript annotations from UniProt, etc"
    def __init__(self, annotSet, transId):
        self.transId = transId

class TargetTranscript:
    def __init__(self, transId):
        self.transId = transId

def _intersectAnnots(srcAnnotMappings, targetAnnotMappings):
    pass

def _intersectTransAnnots(srcTransAnnotMappings, targetTransAnnotMappings):
    pass

def _annotCmpTransUniprot(srcTrans, targetTrans, featurePairings):
    pass
