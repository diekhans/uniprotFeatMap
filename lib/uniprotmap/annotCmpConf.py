"""
Configuration objects for cross-species annotation comparisons.
These are loaded from a python file.
"""
from dataclasses import dataclass
from typing import List
from pycbio.sys.symEnum import SymEnum, auto

class AnnotSet(SymEnum):
    """Set of annotations"""
    UniProt = auto()
    TrEMBL = auto()
    UniProtMapped = auto()
    InterPro = auto()

@dataclass
class SrcAnnotConf:
    """Source (reference) annotation data to use in comparison, mapped to the target assembly"""
    annotSet: AnnotSet
    uniprotAnnotTsv: str
    annot2GenomePslFile: str
    annot2GenomeRefTsv: str
    xspeciesTrans2TransPslFile: str

@dataclass
class TargetAnnotConf:
    """Target annotation to compare with source"""
    annotSet: AnnotSet
    annotTsv: str    # either uniprotAnnotTsv or interproAnnotTsv output
    annot2TransPsl: str
    annot2GenomePsl: str
    annot2GenomeRefTsv: str

@dataclass
class AnnotCmpConf:
    """Configuration of comparisons"""
    srcAnnotConf: SrcAnnotConf
    targetAnnotConfs: List[TargetAnnotConf]
