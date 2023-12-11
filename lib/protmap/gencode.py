from collections import defaultdict
from pycbio.tsv import TsvReader
from pycbio.hgdata.psl import Psl, PslBlock, PslReader
from pycbio.hgdata.genePred import GenePredReader
from pycbio.hgdata.rangeFinder import RangeFinder

# Note: multiple annotation of a transcript is not needed in current
# GenCode, but will be needed for RefSeq, so leave logic here.

def isTargetedTranscriptType(gencodeMeta):
    "is this transcript one that we should map with swissport?"
    return gencodeMeta.transcriptType == "protein_coding"

class GencodeMetaTbl(list):
    def __init__(self, gencodeMetaTsv):
        self.byTranscriptId = {}
        self.codingTransIds = set()
        for row in TsvReader(gencodeMetaTsv):
            self._readRow(row)

    def _readRow(self, row):
        self.byTranscriptId[row.transcriptId] = row
        if isTargetedTranscriptType(row):
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
        self._loadAnnots(gencodeGp)
        self.byTransId.default_factory = None
        self.rangeIdx = None

    def findEntries(self, transId):
        return self.byTransId.get(transId, ())

    def getEntries(self, transId, chrom):
        entries = self.findEntries(transId, chrom)
        if len(entries) is None:
            raise Exception(f"GENCODE not found for '{transId}'")
        return entries

    def findEntry(self, transId, chrom):
        # handle multi-location entries
        entries = self.byTransId.get(transId, ())
        for entry in entries:
            if entries.alignPsl.tName == chrom:
                return entry
        return None

    def getEntry(self, transId, chrom):
        entry = self.findEntry(transId, chrom)
        if entry is None:
            raise Exception(f"GENCODE not found for '{transId}' on chrom '{chrom}'")
        return entry

    def _loadAligns(self, gencodePsl):
        for psl in PslReader(gencodePsl):
            self.byTransId[psl.qName].append(GencodeDataTbl._Entry(psl))

    def _loadAnnots(self, gencodeGp):
        for gp in GenePredReader(gencodeGp):
            entry = self.getEntry(gencodeGp.name, gencodeGp.chrom)
            entry.annotGp = gp

    def _buildRangeIdx(self):
        self.rangeIdx = RangeFinder()
        for entries in self.byTransId.values():
            for entry in entries:
                self.rangeIdx.add(entry.alignPsl.tName, entry.alignPsl.tStart, entry.alignPsl.tEnd)

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


##
# build PSL with UTR unaligned
##
def _makeCdsPslBlock(gp, blk):
    tStart = max(gp.cdsStart, blk.tStart)
    size = min(gp.cdsEnd, blk.tEnd) - tStart
    if size <= 0:
        return None
    qStart = blk.qStart + (tStart - blk.tStart)
    return PslBlock(qStart, tStart, size)

def _addCdsPslBlocks(gp, psl, cdsPsl):
    for blk in psl.blocks:
        cdsBlk = _makeCdsPslBlock(gp, blk)
        if cdsBlk is not None:
            cdsPsl.addBlock(cdsBlk)

def gencodeMakeCdsPsl(gp, psl):
    "modify a transcript PSL to make the UTR unaligned, returning a new PSL"
    cdsPsl = Psl.create(qName=psl.qName, qSize=psl.qSize,
                        tName=psl.tName, tSize=psl.tSize,
                        strand=psl.strand)
    _addCdsPslBlocks(gp, psl, cdsPsl)

    # fixup bounds
    cdsPsl.tStart = cdsPsl.blocks[0].tStartPlus
    cdsPsl.tEnd = cdsPsl.blocks[-1].tEndPlus
    cdsPsl.qStart = cdsPsl.blocks[0].qStartPlus
    cdsPsl.qEnd = cdsPsl.blocks[-1].qEndPlus
    return cdsPsl
