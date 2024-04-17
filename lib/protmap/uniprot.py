"""
Reads files create by uniprotToTab, and other uniport support
"""

from collections import defaultdict
from pycbio.sys.symEnum import SymEnum, auto
from pycbio.tsv import TsvReader, TsvRow
from protmap import dropVersion

# WARNING: UniProt is 1-based, open-end

class UniProtDataSet(SymEnum):
    SwissProt = 1
    TrEMBL = 2

class UniProtCategory(SymEnum):
    "categories assign my Max in genomic UniProt tracks"
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

class TransMatchStatus(SymEnum):
    same_version = auto()
    diff_version = auto()
    other_isoform = auto()

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
        raise KeyError(f"key '{key}' already in unique index")
    idx[key] = value

def addValuesToMultiIdx(idx, keys, value):
    for key in keys:
        idx[key].append(value)

class UniProtMeta(TsvRow):
    "one row record of UniProt metadata"

    def __init__(self, reader, row):
        super().__init__(reader, row)
        # split comma list fields
        self.ensemblGeneIds = splitMetaList(self.ensemblGene)
        self.ensemblGeneAccs = splitDropVersion(self.ensemblGene)
        self.ensemblTransIds = splitMetaList(self.ensemblTrans)
        self.ensemblTransAccs = splitDropVersion(self.ensemblTrans)
        self.uniprotIsoIds = splitMetaList(self.isoIds)

    def isCanonProtTrans(self, transId):
        "does this UniProt reference this transcript?"
        return dropVersion(transId) in self.ensemblTransAccs


class UniProtMetaTbl:
    """reads swissprot.9606.tab, trembl.9606.tab"""
    def __init__(self, uniprotMetaTsv):
        self.byAcc = {}
        self.byMainIsoAcc = {}
        self.byGeneName = defaultdict(list)
        self.byGeneAcc = defaultdict(list)
        self.byTranscriptAcc = defaultdict(list)
        self.byUniProtIsoId = defaultdict(list)
        for row in TsvReader(uniprotMetaTsv, rowClass=UniProtMeta):
            self._readRow(row)
        self.byGeneName.default_factory = None
        self.byGeneAcc.default_factory = None
        self.byTranscriptAcc.default_factory = None
        self.byUniProtIsoId.default_factory = None

    def _readRow(self, row):
        addUniqueToIdx(self.byAcc, row.acc, row)
        addUniqueToIdx(self.byMainIsoAcc, row.mainIsoAcc, row)
        self.byGeneName[row.geneName].append(row)
        addValuesToMultiIdx(self.byGeneAcc, row.ensemblGeneAccs, row)
        addValuesToMultiIdx(self.byTranscriptAcc, row.ensemblTransAccs, row)
        addValuesToMultiIdx(self.byUniProtIsoId, row.uniprotIsoIds, row)

    def getByAcc(self, acc):
        "Error if not found"
        try:
            return self.byAcc[acc]
        except KeyError as ex:
            raise KeyError(f"UniProt acc not found {acc}") from ex

    def getByIsoId(self, isoId):
        "Error if not found"
        try:
            return self.byIsoId[isoId]
        except KeyError as ex:
            raise KeyError(f"UniProt isoId not found {isoId}") from ex

    def getTransMeta(self, transAcc):
        "None if not found"
        protMetas = self.byTranscriptAcc.get(transAcc)
        if protMetas is None:
            return None
        if len(protMetas) > 1:
            raise Exception(f"more than one UniProt protein for transcript {transAcc}, found: {protMetas.mainIsoAcc}")
        return protMetas[0]

    def getGeneAccMetas(self, geneAcc):
        "list or None if not found"
        return self.byGeneAcc.get(geneAcc)

    def getGeneNameMetas(self, geneName):
        "list or None if not found"
        return self.byGeneName.get(geneName)

    def isCanonProtTrans(uniprotMetaTbl, transId):
        "does any UniProt reference this transcript?"
        return dropVersion(transId) in uniprotMetaTbl.byTranscriptAcc

    def getCanonEnsemblAccSet(self):
        return frozenset([dropVersion(transId) for transId in self.byTranscriptAcc.keys()])

def uniprotParseAnnotId(annotId):
    "parse the annotation id into its parts: 'Q9BXI3#0' -> ('Q9BXI3' 0)"
    parts = annotId.rsplit('#', 1)
    if len(parts) != 2:
        raise Exception(f"invalid annotation id: '{annotId}'")
    return parts

def uniprotAnnotIdToAcc(annotId):
    "parse the annotation id into accession: 'Q9BXI3#0' -> 'Q9BXI3'"
    return uniprotParseAnnotId(annotId)[0]

class UniProtAnnotTbl(list):
    """reads swissprot.9606.annots.tab, trembl.9606.annots.tab"""
    def __init__(self, uniprotAnnotsTsv):
        self.byAnnotId = {}
        mainIsoNextId = defaultdict(int)
        for row in TsvReader(uniprotAnnotsTsv, typeMap={"begin": int, "end": int}):
            self._readRow(row, mainIsoNextId)

    def _readRow(self, row, mainIsoNextId):
        # Add a unique, reproducible annotation id in the form mainIsoAcc#annotIdx, where the annotIdx is
        # the relative index of the annotation for that acc.
        annotId = row.mainIsoAcc + '#' + str(mainIsoNextId[row.mainIsoAcc])
        mainIsoNextId[row.mainIsoAcc] += 1

        setattr(row, "annotId", annotId)
        self.byAnnotId[annotId] = row

        # FIXME: fix feature type due to parser not set featureType to "phosphorylation site"
        # from "modified residue" as it does other types
        if row.shortFeatType == "phos":
            row.featType = "phosphorylation site"

        self.append(row)

    def getByAnnotId(self, annotId):
        annot = self.byAnnotId.get(annotId)
        if annot is None:
            raise Exception(f"annotId '{annotId}' not found in annotation table")
        return annot
