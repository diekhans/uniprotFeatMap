"""
Abstraction of a gene set, allowing code to work with GENCODE, RefSeq and CAT.
"""
from collections import namedtuple, defaultdict
from pycbio.sys.symEnum import SymEnum, auto
from pycbio.hgdata.psl import PslReader
from pycbio.hgdata.genePred import GenePredReader
from pycbio.hgdata.rangeFinder import RangeFinder
from uniprotmap import dropVersion

##
# Note: this came from a different project, so not all of the functionality is
# used.
#
# If metadata is not need, the GeneSetData class and load functions maybe used
# directly.
##

class GeneSetName(SymEnum):
    """Identifies the supported gene sets"""
    GENCODE = auto()
    CAT1 = auto()

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

    def haveTranscriptAcc(self, transAcc):
        return transAcc in self.byTransAcc

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

    def getGeneByTranscriptId(self, transId):
        return self.getTranscript(transId).gene

    def getGeneTranscripts(self, geneId):
        try:
            return self.byGeneId[geneId]
        except Exception as ex:
            raise GeneSetError(f"can't find gene for gene `{geneId}'") from ex

class Entry:
    "PSL and option genePred"
    __slots__ = ("psl", "gp")

    def __init__(self):
        self.psl = None
        self.gp = None

    @property
    def name(self):
        return self.psl.qName if self.psl is not None else self.gp.name

    @property
    def chrom(self):
        return self.psl.tName if self.psl is not None else self.gp.chrom

    @property
    def start(self):
        return self.psl.tStart if self.psl is not None else self.gp.txStart

    @property
    def end(self):
        return self.psl.tEnd if self.psl is not None else self.gp.txEnd

    @property
    def strand(self):
        return self.psl.qStrand if self.psl is not None else self.gp.strand

class GeneSetData:
    """contains PSL alignments and optional genePred annotation records."""

    ##
    # Note: RefSeq and older GENCODEs will have same transcript ids for PARs, so we must keep
    # list of entries
    ##

    def __init__(self):
        self.entries = []
        self.byTransId = defaultdict(list)
        self.rangeIdx = None

    def finish(self):
        self.byTransId.default_factory = None

    def _obtainEntry(self, name, chrom):
        entry = self.findEntry(name, chrom)
        if entry is None:
            entry = Entry()
            self.entries.append(entry)
            self.byTransId[name].append(entry)
        return entry

    def addAlign(self, psl):
        self._obtainEntry(psl.qName, psl.tName).psl = psl

    def addAnnot(self, gp):
        self._obtainEntry(gp.name, gp.chrom).gp = gp

    def _buildRangeIdx(self):
        self.rangeIdx = RangeFinder()
        for entries in self.byTransId.values():
            for entry in entries:
                self.rangeIdx.add(entry.name, entry.start, entry.end)

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
            if entry.chrom == chrom:
                return entry
        return None

    def getEntry(self, transId, chrom):
        entry = self.findEntry(transId, chrom)
        if entry is None:
            raise GeneSetError(f"transcript not found for '{transId}' on chrom '{chrom}'")
        return entry

    def getAlign(self, transId, chrom):
        return self.getEntry(transId, chrom).psl

    def getAnnot(self, transId, chrom):
        return self.getEntry(transId, chrom).gp

    def getOverEntries(self, name, start, end, strand):
        if self.rangeIdx is None:
            self._buildRangeIdx()
        for entry in self.rangeIdx.overlapping(self, name, start, end):
            if entry.strand == strand:
                yield entry

def geneSetLoadAnnotPsl(geneSetData, trans2GenomePslFile):
    """"load PSLs into a GeneSetData object"""
    for psl in PslReader(trans2GenomePslFile):
        geneSetData.addAlign(psl)

def geneSetLoadAnnotGp(geneSetData, annotGpFile):
    """"load genePreds into a GeneSetData object"""
    for gp in GenePredReader(annotGpFile):
        geneSetData.addAnnot(gp)


class GeneSet:
    """Meta data and data for GeneSet. A factory function builds this
    to allow set specific functions. Data that is loaded depends on
    what is requested by the program."""
    def __init__(self, geneSetName):
        self.geneSetName = geneSetName
        self.meta = GeneSetMetadata()
        self.data = GeneSetData()
        self.transFa = None

    def finish(self):
        self.data.finish()
        self.meta.finish()


def geneSetFactory(geneSetName, *, geneSetMetadata=None, trans2GenomePslFile=None, transGenomeGpFile=None, transFa=None):
    """Build gene set object of the specified type GeneSet"""

    if geneSetName is GeneSetName.GENCODE:
        from uniprotmap.gencode import gencodeGeneSetFactory
        return gencodeGeneSetFactory(geneSetName, geneSetMetadata=geneSetMetadata,
                                     trans2GenomePslFile=trans2GenomePslFile, transGenomeGpFile=transGenomeGpFile,
                                     transFa=transFa)
    elif geneSetName is GeneSetName.CAT1:
        from uniprotmap.catgenes import cat1GeneSetFactory
        return cat1GeneSetFactory(geneSetName, geneSetMetadata=geneSetMetadata,
                                  trans2GenomePslFile=trans2GenomePslFile, transGenomeGpFile=transGenomeGpFile, transFa=transFa)
    else:
        assert False, f"Bug: no handler for {geneSetName}"
