"""
Metadata support associated with annotation mappings.

annotMapIds are in the form: <canon_acc>|<annot_idx>|<map_idx>
"""

from collections import defaultdict
from pycbio.sys import fileOps
from pycbio.hgdata.coords import Coords
from pycbio.tsv import TsvReader, TsvRow, strOrNoneType, intOrNoneType
from uniprotmap import annotMapIdFmt, annotMapIdToAnnotId

###

class AnnotTransRef(TsvRow):
    """ Metadata for mapping of an annotation to a transcript.
    This parallels the PSL file containing the annotation to transcript
    alignments.  Annotations that did not map to the transcript have entries
    that will indicate their relative position withing the protein but no
    alignment index.

    Attributes:
        annotId: identifies annotation within an UniProt entry (Q9BXI3|0) [not stored in file]
        annotMapId: identifies mapping of annotation to a given transcript (Q9BXI3|0|0), None if
            not mapped.
        transcriptPos: chr1:11674479-11691650 position of transcript (FIXME needed, in PSL)
        transcriptId: ENST00000235310.7
        xspeciesSrcTransId: transcript id of the source transcript that was mapped to the other
            species transcript
        alignIdx: zero-based line numner of PSL alignment file associated with this annotation, or
           None if not mapped
    """
    __slots__ = ()

    def __init__(self, reader, row):
        super().__init__(reader, row)
        self.annotId = annotMapIdToAnnotId(self.annotMapId)

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


class AnnotTransRefWriter:
    header = ["annotMapId", "transcriptPos", "transcriptId", "xspeciesSrcTransId", "alignIdx"]

    def __init__(self, annotTransRefTsv):
        self.fh = fileOps.opengz(annotTransRefTsv, 'w')
        self.idxCounter = defaultdict(int)  # used to get globally unique id
        fileOps.prRow(self.fh, self.header)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    def close(self):
        self.fh.close()
        self.fh = None

    def __del__(self):
        if self.fh is not None:
            self.close()

    def write(self, annotId, transcriptId, transcriptPos, xspeciesSrcTransId, alignIdx):
        "write record, return assigned annotMapId"
        annotMapId = annotMapIdFmt(annotId, self.idxCounter[annotId])
        self.idxCounter[annotId] += 1
        fileOps.prRowv(self.fh, annotMapId, transcriptPos, transcriptId, alignIdx, xspeciesSrcTransId)


def xrefToItemArgs(annotTransRef):
    "convert xref into into [name, start, end]"
    coords = Coords.parse(annotTransRef.transcriptPos)
    return annotTransRef.transcriptId, coords.start, coords.end
