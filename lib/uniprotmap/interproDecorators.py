"""
Support for creating InterProScan decorators
"""

from pycbio.sys.svgcolors import SvgColors
from pycbio.sys.color import Color
from pycbio.sys.symEnum import SymEnum, auto
from uniprotmap import dropVersion
from pycbio.hgdata.bed import encodeRow
from pycbio.hgdata.decoration import Decoration

###
# Color logic.  This mostly matches Interpro genomic track, however colors for
# overlay are not visible with track colors, so these are adjusted.
###

def _mkcolor(r, g, b, a=None):
    return Color.fromRgb8(r, g, b, a)

INTERPRO_COLOR = _mkcolor(0, 191, 255)    # deepskyblue
OTHER_COLOR = _mkcolor(100, 100, 0),      # olive

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
