"""
File that cross-references mapped annotations to transcript metadata
"""

import os
from pycbio.sys import fileOps
os.environ["OPENBLAS_NUM_THREADS"] = "1"
import pandas as pd

class AnnotTransRefs:
    """look up transcript for a PSL row"""
    def __init__(self, annotTransRefTsv):
        self.df = pd.read_table(annotTransRefTsv)
        self.df.set_index('alignIdx', inplace=True, drop=False, verify_integrity=True)

    def get(self, annotId, alignIdx):
        annotTransRef = self.df[self.df.index == alignIdx].iloc[0]
        if annotId != annotTransRef.annotId:
            raise Exception(f"annotGenomePsl '{annotId}' and  annotTransRefTsv '{annotTransRef.annotId}' out-of-sync")
        return annotTransRef

def annotTransRefCreate(annotTransRefTsv):
    "create a new ref TSV, and write header"
    annotTransRefFh = fileOps.opengz(annotTransRefTsv, 'w')
    fileOps.prRowv(annotTransRefFh, "annotId", "transcriptPos", "transcriptId", "alignIdx")
    return annotTransRefFh

def annotTransRefWrite(annotTransRefFh, annotId, transcriptPos, transcriptId, alignIdx):
    fileOps.prRowv(annotTransRefFh, annotId, transcriptPos, transcriptId, alignIdx)
