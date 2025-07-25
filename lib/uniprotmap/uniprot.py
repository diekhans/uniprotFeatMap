"""
Reads files create by uniprotToTab, and other uniport support
"""

from collections import defaultdict
from pycbio.sys.symEnum import SymEnum, auto
from pycbio.tsv import TsvReader, TsvRow
from uniprotmap import dropVersion, annotIdFmt

# WARNING: UniProt is 1-based, open-end

class UniProtError(Exception):
    pass

class UniProtDataSet(SymEnum):
    SwissProt = 1
    TrEMBL = 2

class UniProtCategory(SymEnum):
    "categories assign by Max in genomic UniProt tracks"
    Mut = auto()
    Splice = auto()
    Struct = auto()
    LocTransMemb = auto()
    LocExtra = auto()
    LocCytopl = auto()
    LocSignal = auto()
    Chain = auto()
    Repeat = auto()
    Conflict = auto()
    DisulfBond = auto()
    Domain = auto()
    Modif = auto()
    Interest = auto()
    Other = auto()

class TransCategory(SymEnum):
    canonical = auto()
    canonical_diff_version = auto()
    noncanonical = auto()

def splitMetaList(val):
    """ split strings like: ENST00000369413.8|ENST00000528909.1"""
    if val == "":
        return []
    return val.split('|')

def splitDropVersion(idents):
    return [dropVersion(id) for id in splitMetaList(idents)]

def dropUniportIsoformModifier(acc):
    return acc.split('-')[0]

def addUniqueToIdx(idx, key, value):
    if key in idx:
        raise UniProtError(f"key '{key}' already in unique index")
    idx[key] = value

def addValuesToMultiIdx(idx, keys, value):
    for key in keys:
        idx[key].append(value)

class UniProtMeta(TsvRow):
    """one row record of UniProt metadata

    acc dataset mainIsoAcc orgName orgCommon taxonId name accList
    protFullNames protShortNames protAltFullNames protAltShortNames geneName
    geneSynonyms isoNames geneOrdLocus geneOrf hgncSym hgncId refSeq
    refSeqProt entrezGene ensemblGene ensemblProt ensemblTrans kegg emblMrna
    emblMrnaProt emblDna emblDnaProt pdb ec uniGene omimGene omimPhenotype
    subCellLoc functionText isoIds
    """

    def __init__(self, reader, row):
        super().__init__(reader, row)
        # split '|' separated list fields
        self.geneNames = frozenset(splitMetaList(self.geneName))
        self.ensemblGeneIds = frozenset(splitMetaList(self.ensemblGene))
        self.ensemblGeneAccs = frozenset(splitDropVersion(self.ensemblGene))
        self.ensemblTransIds = frozenset(splitMetaList(self.ensemblTrans))
        self.ensemblTransAccs = frozenset(splitDropVersion(self.ensemblTrans))
        self.otherIsoIds = frozenset(splitMetaList(self.isoIds))

    def isCanonProtTrans(self, transId):
        "does this UniProt reference this transcript?"
        return dropVersion(transId) in self.ensemblTransAccs


class UniProtMetaTbl(list):
    """reads swissprot.9606.tab, trembl.9606.tab"""
    def __init__(self, uniprotMetaTsv):
        self.byAcc = {}
        self.byMainIsoAcc = {}
        self.byGeneName = defaultdict(list)
        self.byGeneAcc = defaultdict(list)
        self.byTranscriptAcc = defaultdict(list)
        for row in TsvReader(uniprotMetaTsv, rowClass=UniProtMeta):
            self._readRow(row)
        self.byGeneName.default_factory = None
        self.byGeneAcc.default_factory = None
        self.byTranscriptAcc.default_factory = None

    def _readRow(self, row):
        self.append(row)
        addUniqueToIdx(self.byAcc, row.acc, row)
        addUniqueToIdx(self.byMainIsoAcc, row.mainIsoAcc, row)
        self.byGeneName[row.geneName].append(row)
        addValuesToMultiIdx(self.byGeneAcc, row.ensemblGeneAccs, row)
        addValuesToMultiIdx(self.byTranscriptAcc, row.ensemblTransAccs, row)

    def getByAcc(self, acc):
        "Error if not found"
        try:
            return self.byAcc[acc]
        except Exception as ex:
            raise UniProtError(f"UniProt acc not found {acc}") from ex

    def getByMainIsoAcc(self, isoAcc):
        "Error if not found"
        try:
            return self.byMainIsoAcc[isoAcc]
        except Exception as ex:
            raise UniProtError(f"UniProt main isoform not found {isoAcc}") from ex

    def getTransMeta(self, transAcc):
        "None if not found"
        protMetas = self.byTranscriptAcc.get(transAcc)
        if protMetas is None:
            return None
        if len(protMetas) > 1:
            raise UniProtError(f"more than one UniProt protein for transcript {transAcc}, found: {protMetas.mainIsoAcc}")
        return protMetas[0]

    def getGeneAccMetas(self, geneAcc):
        return self.byGeneAcc.get(geneAcc, ())

    def getGeneNameMetas(self, geneName):
        return self.byGeneName.get(geneName, ())

    def isCanonProtTrans(uniprotMetaTbl, transId):

        "does any UniProt reference this transcript?"
        return dropVersion(transId) in uniprotMetaTbl.byTranscriptAcc

    def getCanonEnsemblAccSet(self):
        return frozenset([dropVersion(transId) for transId in self.byTranscriptAcc.keys()])


class UniprotAnnot(TsvRow):
    """one interpro annotation record"""

    def short(self):
        desc = f"{self.acc}: {self.shortFeatType}"
        if self.comment is not None:
            desc += f"/{self.comment}"
        return desc


class UniProtAnnotTbl(list):
    """reads swissprot.9606.annots.tab or trembl.9606.annots.tab

    acc mainIsoAcc varId featType shortFeatType begin end origAa mutAa dbSnpId
    disRelated disease disCode pmid longName shortName syns subCellLoc comment
    """
    def __init__(self, uniprotAnnotsTsv):
        self.byAnnotId = {}
        self.byMainIsoAcc = defaultdict(list)
        nextFeatId = defaultdict(int)
        for row in TsvReader(uniprotAnnotsTsv, typeMap={"begin": int, "end": int},
                             rowClass=UniprotAnnot):
            self._readRow(row, nextFeatId)

    def _readRow(self, row, nextFeatId):
        annotId = annotIdFmt(row.mainIsoAcc, nextFeatId[row.mainIsoAcc])
        nextFeatId[row.mainIsoAcc] += 1
        row.annotId = annotId
        self.append(row)
        self.byAnnotId[annotId] = row
        self.byMainIsoAcc[row.mainIsoAcc].append(row)

    def getByAnnotId(self, annotId):
        annot = self.byAnnotId.get(annotId)
        if annot is None:
            raise UniProtError(f"UniProt annotId '{annotId}' not found in annotation table")
        return annot
