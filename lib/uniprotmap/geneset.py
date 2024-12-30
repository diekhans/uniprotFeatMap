"""
Abstraction of a gene set, allowing code to work with GENCODE, RefSeq and CAT.
"""
from collections import namedtuple, defaultdict
from pycbio.hgdata.psl import PslReader
from pycbio.hgdata.genePred import GenePredReader
from pycbio.hgdata.rangeFinder import RangeFinder
from uniprotmap import dropVersion

##
# note: this came from a different project, so not all of the functionality is used.
##

class GeneSetError(Exception):
    pass

class GeneMetadata(namedtuple("GeneMetadata",
                              ("geneId",
                               "geneAcc",
                               "geneSymbol",
                               "geneType"))):
    """ metadata for a gene. Accession drops version"""
    pass

class TranscriptMetadata(namedtuple("TranscriptMetadata",
                                    ("transId",
                                     "transAcc",
                                     "transType",
                                     "proteinId",
                                     "gene"))):
    """Metadata for a transcript. The transType is specific to
    the derived type.  Accession drops version.  proteinId can
    be None."""
    pass


class GeneSetMetadata:
    """metadata for a gene set. Only protein coding transcripts should
    be included"""
    def __init__(self):
        self.genes = []
        self.byGeneId = {}
        self.byGeneAcc = {}
        self.byGeneSymbol = {}
        self.transcripts = []
        self.byTransId = {}
        self.byTransAcc = {}
        self.byProteinId = {}
        self.transesByGeneId = defaultdict(list)

        self.geneIds = None
        self.geneAccs = None
        self.transcriptIds = None
        self.transcriptAccs = None

    def finish(self):
        "call after loaded"
        self.transesByGeneId.default_factory = None

        self.geneIds = frozenset(self.byGeneId.keys())
        self.geneAccs = frozenset(self.byGeneAcc.keys())
        self.transcriptIds = frozenset(self.byTransId.keys())
        self.transcriptAccs = frozenset(self.byTransAcc.keys())

    def _addGene(self, geneId, geneSymbol, geneType):
        geneMd = GeneMetadata(geneId, dropVersion(geneId), geneSymbol, geneType)
        self.genes.append(geneMd)
        self.byGeneId[geneMd.geneId] = geneMd
        self.byGeneAcc[geneMd.geneAcc] = geneMd
        self.byGeneSymbol[geneMd.geneSymbol] = geneMd
        return geneMd

    def _obtainGene(self, geneId, geneSymbol, geneType):
        geneMd = self.byGeneId.get(geneId)
        if geneMd is None:
            geneMd = self._addGene(geneId, geneSymbol, geneType)
        return geneMd

    def addTranscript(self, geneId, geneSymbol, geneType, transId, transType, proteinId):
        if transId in self.byTransId:
            raise GeneSetError(f"duplicate transcripts id `{transId}'")
        if proteinId in self.byProteinId:
            raise GeneSetError(f"duplicate protein id `{proteinId}'")

        geneMd = self._obtainGene(geneId, geneSymbol, geneType)
        transMd = TranscriptMetadata(transId, dropVersion(transId), transType,
                                     proteinId, geneMd)
        self.transcripts.append(transMd)
        self.byTransId[transMd.transId] = transMd
        self.byTransAcc[transMd.transAcc] = transMd
        if proteinId is not None:
            self.byProteinId[proteinId] = transMd
        self.transesByGeneId[geneId].append(transMd)
        return transMd

    def getTranscript(self, transId):
        try:
            return self.byTransId[transId]
        except Exception as ex:
            raise GeneSetError(f"can't find metadata for transcript `{transId}'") from ex

    def getTranscriptByAcc(self, transAcc):
        try:
            return self.byTransAcc[transAcc]
        except Exception as ex:
            raise GeneSetError(f"can't find metadata for transcript accession `{transAcc}'") from ex

    def transcriptIdIter(self):
        "iterator for transcript ids"
        return self.byTransId.keys()

    def transcriptAccIter(self):
        "iterator for transcript accessions"
        return self.byTransAcc.keys()

    def getGene(self, geneId):
        try:
            return self.byGeneId[geneId]
        except Exception as ex:
            raise GeneSetError(f"can't find metadata for gene `{geneId}'") from ex

    def getGeneTranscripts(self, geneId):
        try:
            return self.byGeneId[geneId]
        except Exception as ex:
            raise GeneSetError(f"can't find gene for gene `{geneId}'") from ex

class Entry:
    "PSL and option genePred"
    __slots__ = ("alignPsl", "annotGp")

    def __init__(self, alignPsl):
        self.alignPsl = alignPsl
        self.annotGp = None  # load second

class GeneSetData:
    """contains PSL alignments and optional genePred annotation records."""

    # Note: older GENCODEs and RefSeq will have same transcript ids for PARs, so we must keep
    # list of entries
    ##

    def __init__(self):
        self.entries = []
        self.byTransId = defaultdict(list)
        self.rangeIdx = None

    def finish(self):
        self.byTransId.default_factory = None

    def addAlign(self, psl):
        self.byTransId[psl.qName].append(Entry(psl))

    def addAnnot(self, gp):
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
        if len(entries) == 0:
            raise GeneSetError(f"transcript not found for '{transId}'")
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
            raise GeneSetError(f"transcript not found for '{transId}' on chrom '{chrom}'")
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

def geneSetDataLoad(alignPslFile, annotGenePredFile=None):
    """"load PSL and options genePreds into a GeneSetData object"""
    geneSetData = GeneSetData()

    for psl in PslReader(alignPslFile):
        geneSetData.addAlign(psl)

    if annotGenePredFile is not None:
        for gp in GenePredReader(annotGenePredFile):
            geneSetData.addAnnot(gp)

    geneSetData.finish()
    return geneSetData
