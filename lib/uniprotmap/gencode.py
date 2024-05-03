import sys
import os.path as osp
from collections import defaultdict
from pycbio.tsv import TsvReader
from pycbio.hgdata.psl import PslReader
from pycbio.hgdata.genePred import GenePredReader
from pycbio.hgdata.rangeFinder import RangeFinder

sys.path.insert(0, osp.normpath(osp.join(osp.dirname(__file__), "../lib")))
from uniprotmap import dropVersion

# Note: multiple annotation of a transcript is not needed in current
# GenCode, but will be needed for RefSeq, so leave logic here.

_codingTranscriptTypes = frozenset(["protein_coding",
                                    "nonsense_mediated_decay",
                                    "non_stop_decay",
                                    "protein_coding_LoF",
                                    "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene",
                                    "IG_C_gene", "IG_J_gene", "IG_D_gene", "IG_V_gene"])

def isCodingTranscriptType(gencodeMeta):
    "is this transcript one that we should map with swissport?"
    return gencodeMeta.transcriptType in _codingTranscriptTypes

class GencodeMetaTbl(list):
    def __init__(self, gencodeMetaTsv):
        self.byTranscriptId = {}
        self.byTranscriptAcc = {}
        self.codingTransIds = set()
        for row in TsvReader(gencodeMetaTsv):
            self._readRow(row)

    def _readRow(self, row):
        self.byTranscriptId[row.transcriptId] = row
        self.byTranscriptAcc[dropVersion(row.transcriptId)] = row
        if isCodingTranscriptType(row):
            self.codingTransIds.add(row.transcriptId)

    def _doGetTrans(self, transId):
        # handle PAR_Y ids (e.g. ENST00000359512.8_PAR_Y)
        if transId.endswith("_PAR_Y"):
            transId = transId.split('_', 1)[0]
        try:
            return self.byTranscriptId.get(transId)
        except Exception as ex:
            raise Exception(f"can't find GENCODE metadata for '{transId}'") from ex

    def getTrans(self, transId):
        try:
            return self._doGetTrans(transId)
        except Exception as ex:
            raise Exception(f"getTrans error for '{transId}'") from ex

    def getGeneId(self, transId):
        return self.getTrans(transId).geneId

    def getGeneName(self, transId):
        return self.getTrans(transId).geneName

    def getCodingTransIds(self):
        return self.codingTransIds

    def getCodingTransAccSet(self):
        return frozenset([dropVersion(transId) for transId in self.getCodingTransIds()])

    def getTransType(self, transId):
        return self.getTrans(transId).transcriptType

class GencodeDataTbl:
    class _Entry:
        __slots__ = ("alignPsl", "annotGp")

        def __init__(self, alignPsl):
            self.alignPsl = alignPsl
            self.annotGp = None  # load second

    "all annotations, optionally with genePreds"
    def __init__(self, gencodePsl, gencodeGp=None):
        # older GENCODEs and RefSeq will have same transcript ids for PARs, so a list is required
        self.byTransId = defaultdict(list)
        self._loadAligns(gencodePsl)
        if gencodeGp is not None:
            self._loadAnnots(gencodeGp)
        self.byTransId.default_factory = None
        self.rangeIdx = None

    def _loadAligns(self, gencodePsl):
        for psl in PslReader(gencodePsl):
            self.byTransId[psl.qName].append(GencodeDataTbl._Entry(psl))

    def _loadAnnots(self, gencodeGp):
        for gp in GenePredReader(gencodeGp):
            entry = self.getEntry(gp.name, gp.chrom)
            entry.annotGp = gp

    def _buildRangeIdx(self):
        self.rangeIdx = RangeFinder()
        for entries in self.byTransId.values():
            for entry in entries:
                self.rangeIdx.add(entry.alignPsl.tName, entry.alignPsl.tStart, entry.alignPsl.tEnd)

    def findEntries(self, transId):
        return self.byTransId.get(transId, ())

    def getEntries(self, transId):
        entries = self.findEntries(transId)
        if len(entries) is None:
            raise Exception(f"GENCODE metadata not found for '{transId}'")
        return entries

    def findEntry(self, transId, chrom):
        # handle multi-location entries
        entries = self.byTransId.get(transId, ())
        for entry in entries:
            if entry.alignPsl.tName == chrom:
                return entry
        return None

    def getEntry(self, transId, chrom):
        entry = self.findEntry(transId, chrom)
        if entry is None:
            raise Exception(f"GENCODE metadata not found for '{transId}' on chrom '{chrom}'")
        return entry

    def getAlign(self, transId, chrom):
        return self.getEntry(transId, chrom).alignPsl

    def getAnnot(self, transId, chrom):
        return self.getEntry(transId, chrom).annotGp

    def getOverEntries(self, tName, tStart, tEnd, qStrand):
        if self.rangeIdx is None:
            self._buildRangeIdx()
        for entry in self.rangeIdx.overlapping(self, tName, tStart, tEnd):
            if entry.annotPsl.qStrand == qStrand:
                yield entry
