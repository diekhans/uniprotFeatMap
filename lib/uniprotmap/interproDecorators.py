"""
Support for creating InterProScan decorators
"""

from pycbio.sys.svgcolors import SvgColors
from pycbio.sys.symEnum import SymEnum, auto
from pycbio.hgdata.bed import encodeRow
from pycbio.hgdata.decoration import Decoration

###
# Color logic.  This mostly matches Interpro genomic track, however colors for
# overlay are not visible with track colors, so these are adjusted.
###

INTERPRO_COLOR = SvgColors.deepskyblue
OTHER_COLOR = SvgColors.olive

class AnnotCategory(SymEnum):
    "indicates type of anaylsis"
    InterPro = auto()
    other = auto()

class InterproDecoration(Decoration):
    """Holds a decoration BED + extra data.  This matches etc/interproDecoration.as. """

    # Keep in schema order:
    __slots__ = ("proteinAcc", "analysis", "accession", "description", "annotCategory",
                 "interproAcc", "interproDesc", "signatureAcc", "signatureDesc",)

    def __init__(self, chrom, bedBlocks, name, strand, color, itemName,
                 itemStart, itemEnd, glyph, fillColor, *, proteinAcc,
                 analysis, accession, description, annotCategory, interproAcc,
                 interproDesc, signatureAcc, signatureDesc):
        super(InterproDecoration, self).__init__(chrom, bedBlocks[0].start, bedBlocks[-1].end, name,
                                                 itemName, itemStart, itemEnd,
                                                 strand=strand, itemRgb=color, blocks=bedBlocks,
                                                 glyph=glyph, fillColor=fillColor)
        self.proteinAcc = proteinAcc
        self.analysis = analysis
        self.accession = accession
        self.description = description
        self.annotCategory = annotCategory
        self.interproAcc = interproAcc
        self.interproDesc = interproDesc
        self.signatureAcc = signatureAcc
        self.signatureDesc = signatureDesc

    def toRow(self):
        return super().toRow() + encodeRow([getattr(self, f) for f in self.__slots__])
