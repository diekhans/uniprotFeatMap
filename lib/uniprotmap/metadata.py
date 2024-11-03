"""
Metadata support associated with annotation mappings.

annotMapIds are in the form: <canon_acc>|<annot_idx>|<map_idx>
"""

from pycbio.sys import fileOps
from pycbio.tsv import TsvReader, TsvRow, strOrNoneType, intOrNoneType

###

class AnnotTransRef(TsvRow):
    """ Metadata for mapping of an annotation to a transcript.
    This parallels the PSL file containing the annotation to transcript
    alignments.  Annotations that did not map to the transcript have entries
    that will indicate their relative position withing the protein but no
    alignment index.

    Attributes:
        annotId: identifies annotation within an UniProt entry (Q9BXI3|0)
        annotMapId: identifies mapping of annotation to a given transcript (Q9BXI3|0|0), None if
            not mapped.
        transcriptPos: chr1:11674479-11691650 position of transcript (FIXME needed, in PSL)
        transcriptId: ENST00000235310.7
        xspeciesSrcTransId: transcript id of the source transcript that was mapped to the other
            species transcript
        alignIdx: zero-based line numner of PSL alignment file associated with this annotation, or
           None if not mapped
    """
    # just here for documentation, could be a named tuple
    __slots__ = ()
    pass

_annotTransRefTypeMap = {
    "xspeciesSrcTransId": strOrNoneType,
    "alignIdx": intOrNoneType,
}

class AnnotTransRefError(Exception):
    pass

class AnnotTransRefs:
    """look up transcript for a PSL row"""
    def __init__(self, annotTransRefTsv):
        self.byAlignIdx = {}
        for row in TsvReader(annotTransRefTsv, typeMap=_annotTransRefTypeMap, rowClass=AnnotTransRef):
            self._readRow(row)

    def _readRow(self, row):
        self.byAlignIdx[row.alignIdx] = row

    def get(self, annotId, alignIdx):
        annotTransRef = self.byAlignIdx[alignIdx]
        if annotTransRef.annotId != annotId:
            raise AnnotTransRefError(f"annotGenomePsl '{annotId}' and  annotTransRefTsv '{annotTransRef.annotId}' out-of-sync")
        return annotTransRef

def annotTransRefOpen(annotTransRefTsv):
    "create a new ref TSV, and write header"
    annotTransRefFh = fileOps.opengz(annotTransRefTsv, 'w')
    hdr = ["annotId", "annotMapId", "transcriptPos", "transcriptId", "xspeciesSrcTransId", "alignIdx"]
    fileOps.prRow(annotTransRefFh, hdr)
    return annotTransRefFh

def annotTransRefWrite(annotTransRefFh, annotId, annotMapId, transcriptPos, transcriptId, xspeciesSrcTransId, alignIdx):
    row = [annotId, annotMapId, transcriptPos, transcriptId, alignIdx, xspeciesSrcTransId]
    fileOps.prRow(annotTransRefFh, row)
