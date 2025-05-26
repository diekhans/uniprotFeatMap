"""
Metadata support associated with annotation mappings.

annotMapIds are in the form: <canon_acc>|<annot_idx>|<map_idx>
"""

from collections import defaultdict, namedtuple
from pycbio.sys import fileOps
from pycbio.hgdata.coords import Coords
from pycbio.tsv import TsvReader, strOrNoneType, intOrNoneType
from uniprotmap import annotMapIdFmt, annotMapIdToAnnotId

###

_annotTransRefRowFields = ("annotMapId", "transcriptPos", "transcriptId",
                           "xspeciesSrcTransId", "alignIdx")

class AnnotTransRef(namedtuple("AnnotTransRef",
                    ("annotId",) + _annotTransRefRowFields)):
    """ Metadata for mapping of an annotation to a transcript.
    This parallels the PSL file containing the annotation to transcript
    alignments.  Annotations that did not map to the transcript have entries
    that will indicate their relative position withing the protein but no
    alignment index.

    Attributes:
        annotId: identifies annotation within an UniProt entry (Q9BXI3|0) [not stored in file]
        annotMapId: identifies mapping of annotation to a given transcript (Q9BXI3|0|0), None if
            not mapped.
        transcriptPos: position of transcript as a Coords object
        transcriptId: ENST00000235310.7
        xspeciesSrcTransId: transcript id of the source transcript that was mapped to the other
            species transcript
        alignIdx: zero-based line numner of PSL alignment file associated with this annotation, or
           None if not mapped
    """
    # This needs to be a namedtuple to be pickled
    __slots__ = ()

def _annotTransRefParseRow(reader, row):
    values = [row[reader.colMap[f]] for f in _annotTransRefRowFields]
    annotId = annotMapIdToAnnotId(values[0])
    return AnnotTransRef(annotId, *values)

_annotTransRefTypeMap = {
    "xspeciesSrcTransId": strOrNoneType,
    "alignIdx": intOrNoneType,
    "transcriptPos": Coords.parse,
}

class AnnotTransRefError(Exception):
    pass

def annotTransRefReader(annot2TransRefTsv):
    return TsvReader(annot2TransRefTsv, typeMap=_annotTransRefTypeMap, rowClass=_annotTransRefParseRow)

class AnnotTransRefs:
    """look up transcript for a PSL row"""
    def __init__(self, annot2TransRefTsv):
        self.byAlignIdx = {}
        for row in annotTransRefReader(annot2TransRefTsv):
            self._readRow(row)

    def _readRow(self, row):
        self.byAlignIdx[row.alignIdx] = row

    def get(self, annotId, alignIdx):
        annotTransRef = self.byAlignIdx[alignIdx]
        if annotTransRef.annotId != annotId:
            raise AnnotTransRefError(f"annot2GenomePsl '{annotId}' and  annot2TransRefTsv '{annotTransRef.annotId}' out-of-sync")
        return annotTransRef


class AnnotTransRefWriter:
    header = ["annotMapId", "transcriptPos", "transcriptId", "xspeciesSrcTransId", "alignIdx"]

    def __init__(self, annot2TransRefTsv):
        self.fh = fileOps.opengz(annot2TransRefTsv, 'w')
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
        """Write record, return assigned annotMapId. Must be grouped and sorted by mapped transcript id"""
        annotMapId = annotMapIdFmt(annotId, self.idxCounter[annotId])
        self.idxCounter[annotId] += 1
        fileOps.prRowv(self.fh, annotMapId, transcriptPos, transcriptId, alignIdx, xspeciesSrcTransId)

def xrefToItemArgs(annotTransRef):
    "convert xref into into [name, start, end] for decorator"
    return annotTransRef.transcriptId, annotTransRef.transcriptPos.start, annotTransRef.transcriptPos.end
