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

_annot2GenomeRefRowFields = ("annotMapId", "transcriptPos", "transcriptId",
                             "xspeciesSrcTransId", "alignIdx")

class Annot2GenomeRef(namedtuple("Annot2GenomeRef",
                                 ("annotId",) + _annot2GenomeRefRowFields)):
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

def _annot2GenomeRefParseRow(reader, row):
    values = [row[reader.columnSpecs.columnMap[f]] for f in _annot2GenomeRefRowFields]
    annotId = annotMapIdToAnnotId(values[0])
    return Annot2GenomeRef(annotId, *values)

_annot2GenomeRefTypeMap = {
    "xspeciesSrcTransId": strOrNoneType,
    "alignIdx": intOrNoneType,
    "transcriptPos": Coords.parse,
}

class Annot2GenomeRefError(Exception):
    pass

def annot2GenomeRefReader(annot2GenomeRefTsv):
    return TsvReader(annot2GenomeRefTsv, typeMap=_annot2GenomeRefTypeMap, rowClass=_annot2GenomeRefParseRow)

class Annot2GenomeRefs:
    """look up transcript for a PSL row"""
    def __init__(self, annot2GenomeRefTsv):
        self.byAlignIdx = {}
        for row in annot2GenomeRefReader(annot2GenomeRefTsv):
            self._readRow(row)

    def _readRow(self, row):
        self.byAlignIdx[row.alignIdx] = row

    def get(self, annotId, alignIdx):
        annot2GenomeRef = self.byAlignIdx[alignIdx]
        if annot2GenomeRef.annotId != annotId:
            raise Annot2GenomeRefError(f"annot2GenomePsl '{annotId}' and  annot2GenomeRefTsv '{annot2GenomeRef.annotId}' out-of-sync")
        return annot2GenomeRef


class Annot2GenomeRefWriter:
    header = ["annotMapId", "transcriptPos", "transcriptId", "xspeciesSrcTransId", "alignIdx"]

    def __init__(self, annot2GenomeRefTsv):
        self.fh = fileOps.opengz(annot2GenomeRefTsv, 'w')
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

def xrefToItemArgs(annot2GenomeRef):
    "convert xref into into [name, start, end] for decorator"
    return annot2GenomeRef.transcriptId, annot2GenomeRef.transcriptPos.start, annot2GenomeRef.transcriptPos.end
