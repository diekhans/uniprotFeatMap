"""
Reads files create by uniprotToTab, and other uniport support
"""

from collections import defaultdict
from pycbio.sys.svgcolors import SvgColors
from pycbio.sys.symEnum import SymEnum, auto
from pycbio.tsv import TsvReader
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

class UniProtMetaTbl:
    """reads swissprot.9606.tab, trembl.9606.tab"""
    def __init__(self, uniprotMetaTsv):
        self.byAcc = {}
        self.byMainIsoAcc = {}
        self.byGeneName = defaultdict(list)
        self.byGeneAcc = defaultdict(list)
        self.byTranscriptAcc = defaultdict(list)
        self.byIsoId = defaultdict(list)
        for row in TsvReader(uniprotMetaTsv):
            self._readRow(row)
        self.byGeneName.default_factory = None
        self.byGeneAcc.default_factory = None
        self.byTranscriptAcc.default_factory = None
        self.byIsoId.default_factory = None

    def _readRow(self, row):
        setattr(row, 'ensemblGeneIds', splitMetaList(row.ensemblGene))
        setattr(row, 'ensemblGeneAccs', splitDropVersion(row.ensemblGene))
        setattr(row, 'ensemblTransIds', splitMetaList(row.ensemblTrans))
        setattr(row, 'ensemblTransAccs', splitDropVersion(row.ensemblTrans))
        setattr(row, 'isoIds', splitMetaList(row.isoIds))
        addUniqueToIdx(self.byAcc, row.acc, row)
        addUniqueToIdx(self.byMainIsoAcc, row.mainIsoAcc, row)
        self.byGeneName[row.geneName].append(row)
        addValuesToMultiIdx(self.byGeneAcc, row.ensemblGeneAccs, row)
        addValuesToMultiIdx(self.byTranscriptAcc, row.ensemblTransAccs, row)
        addValuesToMultiIdx(self.byIsoId, row.isoIds, row)

    def getByAcc(self, acc):
        "Error if not found"
        try:
            return self.byAcc[acc]
        except KeyError as ex:
            raise KeyError(f"acc not found {acc}") from ex

    def getByIsoId(self, isoId):
        "Error if not found"
        try:
            return self.byIsoId[isoId]
        except KeyError as ex:
            raise KeyError(f"isoId not found {isoId}") from ex

    def getTransMeta(self, transAcc):
        "None if not found"
        protMetas = self.byTranscriptAcc.get(transAcc)
        if protMetas is None:
            return None
        if len(protMetas) > 1:
            raise Exception(f"more than one protein for transcript {transAcc}, found: {protMetas.mainIsoAcc}")
        return protMetas[0]

    def getGeneAccMetas(self, geneAcc):
        "list or None if not found"
        return self.byGeneAcc.get(geneAcc)

    def getGeneNameMetas(self, geneName):
        "list or None if not found"
        return self.byGeneName.get(geneName)

    def getMatchEnsemblAccSet(self):
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
        self.append(row)

    def getByAnnotId(self, annotId):
        annot = self.byAnnotId.get(annotId)
        if annot is None:
            raise Exception(f"annotId '{annotId}' not found in annotation table")
        return annot

######
# The below code is derived from kent/src/hg/utils/otto/uniprot/doUniprot; if
# this converges, we can turn this into a library.
######

def isMutagenesis(annot):
    return annot.featType == "mutagenesis site"

def isVariant(annot):
    return annot.featType == "sequence variant"


###
# color logic to match uniprot track, however colors for overlay are not visible
# with track colors, so these are adjusted. Also use named colors.
###

TREMBLCOLOR = SvgColors.darkseagreen   # was 0,150,250 light blue (dodgerblue)
SWISSPCOLOR = SvgColors.forestgreen    # was 12,12,120 dark blue (navy)

# mapping of annotations columns to colors
featTypeColors = {
    "modified residue": SvgColors.goldenrod,   # was 200,200,0
    "glycosylation site": SvgColors.teal,      # was 0,100,100
    "disulfide bond": SvgColors.dimgray,       # was 100,100,100
    "topological domain": SvgColors.maroon,    # was 100,0,0
    "zinc finger region": SvgColors.olive,     # was 100,100,0
    "transmembrane region": SvgColors.green,   # was 0,150,0
    "signal peptide": SvgColors.deeppink,      # was 255,0,150 light-red
}
commentColor = {
    "Extracellular": SvgColors.darkcyan,       # was 0,110,180 light-blue
    "Cytoplasmic": SvgColors.darkorange,       # was 255,150,0 light orange
}

def getAnnotColorSwissport(annot):
    color = commentColor.get(annot.comment, None)
    if color is None:
        color = featTypeColors.get(annot.featType, None)
    return color if color is not None else SWISSPCOLOR

# TrEMBL categories to provide special coloring
_tremblColorCategories = frozenset((
    UniProtCategory.LocSignal,
    UniProtCategory.LocExtra,
    UniProtCategory.LocTransMemb,
    UniProtCategory.LocCytopl,
))

def getAnnotColorTrembl(annot):
    if getAnnotCategory(annot)(annot)[0] in _tremblColorCategories:
        return getAnnotColorSwissport(annot)
    else:
        return TREMBLCOLOR

def getAnnotColor(annot, dataSet):
    if dataSet == UniProtDataSet.SwissProt:
        return getAnnotColorSwissport(annot)
    else:
        return getAnnotColorTrembl(annot)

def getProblemColor(annot, isTrembl):
    return SvgColors.red

# some feature types should not go into the bed name field.
# for these features, we use the 'comment' as the bed name
# e.g. "region of interest" is not very interesting, the actual
# description is usually much more interesting.
useCommentFeatType = frozenset(["domain", "chain", "region of interest", "topological domain", "short sequence motif"])

# some very common disease codes are not in the diseases.txt file
disShortNames = {
    "hepatocellular carcinoma": "HepatocC",
    "a hepatocellular carcinoma": "HepatocC",
    "a hepatocellular carcinoma sample": "HepatocC",
    "hepatocellular carcinomas": "HepatocC",
    "ovarian cancer": "OvC",
    "colon cancer": "ColC",
    "colorectal cancer": "ColrectC",
    "hepatoblastoma": "HepatoBl",
    "a sporadic cancer": "sporadic",
    "sporadic cancers": "sporadic",
    "non-small cell lung cancer cell lines": "NSCLC-CL",
    "a colon cancer cell line": "ColC-CL"
}

def shortenDisCode(code):
    "shorten some disease names are not shortened yet by UniProt. "
    newCodes = []
    if code != "":
        for splitCode in code.split(" and "):
            newCodes.append(disShortNames.get(splitCode, splitCode))
    return newCodes

def getDiseaseBedName(annot):
    "bed name field for disease mutation"
    disCodes = shortenDisCode(annot.disCode)
    disName = ",".join(disCodes)
    if len(disName) > 30:
        disName = disName[:30] + "..."
    name = "%s-%sdel" % (annot.begin, annot.end)
    if (len(disCodes) > 0) and ("strain" not in disName):
        name += " in %s" % disName
    return name

def getVarationBedName(annot):
    name = "%s%s%s" % (annot.origAa, annot.begin, annot.mutAa)
    if len(name) > 20:
        name = "%daa to %daa" % (len(annot.origAa), len(annot.mutAa))
    return name

def makeChainRangeName(annot):
    "trembl doesn't have comments for chains so make one up"
    # NOT USED FOR DECORATORS
    # doUniprot does this when no comment for chain, however we this doesn't
    # make sense when we have decorators on differ isoforms.  We leave
    # the code here just in case common code is created and getCommentBedName
    # updated to handle both cases
    return "%s:%s-%s" % (annot.acc, str(annot.begin), str(annot.end))

def getCommentBedName(annot):
    name = annot.comment
    if name == "":
        name = annot.featType
    if name == "Intrinsically disordered":
        name = "Disordered"
    return name

def getAnnotDescriptiveName(annot):   # noqa:C901
    "get name to use in BED, based on doUniprot"
    # set the bed name field to a disease, to the mutation or something else
    if (isMutagenesis(annot) or isVariant(annot)):
        if (annot.origAa == ""):
            return getDiseaseBedName(annot)
        else:
            return getVarationBedName(annot)
    if annot.featType in useCommentFeatType:
        return getCommentBedName(annot)
    if annot.featType == "signal peptide":
        return "Signal peptide"
    if annot.featType == "lipid moiety-binding region":
        return "Lipidation"
    if annot.featType == "transmembrane region":
        return "Transmembrane"
    if annot.comment == "Nuclear localization signal":
        return "Nuclear loc"
    if annot.comment == "Involved in receptor recognition and/or post-binding events":
        return "Recept recog"
    if annot.comment == "Fibronectin type-III":
        return "FibronectIII"

    # more general rules
    if annot.comment.startswith("Zinc finger protein "):
        return annot.comment.replace("Zinc finger protein ", "ZF-")
    if annot.comment.startswith("Necessary for the"):
        return annot.comment.replace("Necessary for the ", "")
    if annot.comment.startswith("Interaction with "):
        return annot.comment.replace("Interaction with ", "Int:")
    if annot.comment.startswith("Involved in the"):
        return annot.comment.replace("Involved in the ", "")
    if annot.comment.startswith("Required for "):
        return annot.comment.replace("Required for ", "")
    if annot.comment.startswith("Cleavage; by host "):
        return "Cleave:" + annot.comment.split()[-1]
    if annot.comment == "Receptor-binding motif; binding to human ACE2":
        return "binds ACE2"

    # chain annotations
    if annot.longName != "":
        return annot.longName
    if annot.shortName != "":
        return annot.shortName

    return annot.shortFeatType

def getAnnotShortDescriptiveName(annot):
    "get BED name, adjust if too long"
    name = getAnnotDescriptiveName(annot)
    if len(name) > 17:
        name = name[:14] + "..."
    return name

def getAnnotCategory(annot):  # noqa: C901
    "return a symbolic and description; from uniport trackDb"
    if annot.featType in ["mutagenesis site", "sequence variant"]:
        return (UniProtCategory.Mut, "Mutation")
    elif annot.featType in ["splice variant"]:
        return (UniProtCategory.Splice, "Splice variant")
    elif annot.featType in ["strand", "helix", "coiled-coil region", "turn"]:
        return (UniProtCategory.Struct, "Protein Structure")
    elif annot.featType in ["transmembrane region"]:
        # suggested by Regeneron: create four subtracks, for the subcell. localization indicators
        return (UniProtCategory.LocTransMemb, "Transmembrane Domain")
    elif annot.comment in ["Extracellular"]:
        return (UniProtCategory.LocExtra, "Extracellular Domain")
    elif annot.comment in ["Cytoplasmic"]:
        return (UniProtCategory.LocCytopl, "Cytoplasmic Domains")
    elif annot.featType in ["signal peptide"]:
        return (UniProtCategory.LocSignal, "Signal Peptide")
    elif annot.featType in ["chain"]:
        return (UniProtCategory.Chain, "Polypeptide Chains")
    elif annot.featType in ["repeat"]:
        return (UniProtCategory.Repeat, "Repeat")
    elif annot.featType == "sequence conflict":
        return (UniProtCategory.Conflict, "Sequence Conflict")
    elif annot.featType in ["disulfide bond"]:
        return (UniProtCategory.DisulfBond, "Disulfide Bond")
    elif annot.featType in ["domain", "zinc finger region", "topological domain"]:
        return (UniProtCategory.Domain, "Domain")
    elif annot.featType in ["glycosylation site", "modified residue", "lipid moiety-binding region"]:
        return (UniProtCategory.Modif, "Amino Acid Modification")
    elif annot.featType in ["region of interest"]:
        return (UniProtCategory.Interest, "Region of Interest")
    else:
        return (UniProtCategory.Other, "Other Annotation")

def calcTransMatchStatus(uniprotMeta, transId):
    "determine match of transcript to GENCODE based on what transcripts are listed in metadata"
    if transId in uniprotMeta.ensemblTransIds:
        return TransMatchStatus.same_version
    elif dropVersion(transId) in uniprotMeta.ensemblTransAccs:
        return TransMatchStatus.diff_version
    else:
        return TransMatchStatus.other_isoform

def getColorUses():
    "list of tuples of (color, use_description)"
    uses = [(SWISSPCOLOR, "SwissProt base"),
            (TREMBLCOLOR, "TrEMBL base")]
    for info, color in featTypeColors.items():
        uses.append((color, info))
    for info, color in commentColor.items():
        uses.append((color, info))
    return uses
