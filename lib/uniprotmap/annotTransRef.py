"""
File that cross-references mapped annotations to transcript metadata and alignment position

feature ids are in the form: <uniprot_acc>|<feature_idx>

featAnnotId ids are in the form <uniprot_acc>|<feature_idx>|<annot_idx>`

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
    annotTransRefWrite(annotTransRefFh, "annotId", "transcriptPos", "transcriptId", "alignIdx")
    return annotTransRefFh

def annotTransRefWrite(annotTransRefFh, annotId, transcriptPos, transcriptId, alignIdx):
    fileOps.prRowv(annotTransRefFh, annotId, transcriptPos, transcriptId, alignIdx)
