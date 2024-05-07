"""
File that cross-references mapped annotations to transcript metadata and alignment position

annotIds are in the form: <canon_acc>|<annot_idx>

annotMapIds are in the form: <canon_acc>|<annot_idx>|<map_idx>

"""

from pycbio.sys import fileOps
from pycbio.tsv import TsvReader

class AnnotTransRefError(Exception):
    pass

class AnnotTransRefs:
    """look up transcript for a PSL row"""
    def __init__(self, annotTransRefTsv):
        self.byAlignIdx = {}
        for row in TsvReader(annotTransRefTsv, typeMap={"alignIdx": int}):
            self._readRow(row)

    def _readRow(self, row):
        self.byAlignIdx[row.alignIdx] = row

    def get(self, annotId, alignIdx):
        annotTransRef = self.byAlignIdx[alignIdx]
        if annotTransRef.annotId != annotId:
            raise AnnotTransRefError(f"annotGenomePsl '{annotId}' and  annotTransRefTsv '{annotTransRef.annotId}' out-of-sync")
        return annotTransRef

def annotTransRefCreate(annotTransRefTsv):
    "create a new ref TSV, and write header"
    annotTransRefFh = fileOps.opengz(annotTransRefTsv, 'w')
    fileOps.prRowv(annotTransRefFh, "annotId", "annotMapId", "transcriptPos", "transcriptId", "alignIdx")
    return annotTransRefFh

def annotTransRefWrite(annotTransRefFh, annotId, mapIdx, transcriptPos, transcriptId, alignIdx):
    annotMapId = annotId + '|' + str(mapIdx)
    fileOps.prRowv(annotTransRefFh, annotId, annotMapId, transcriptPos, transcriptId, alignIdx)
