"""
Reads files create by uniprotToTab, and other uniport support
"""

import pandas as pd
from pycbio.sys.svgcolors import SvgColors
from pycbio.sys.symEnum import SymEnum
from protmap import dropVersion, buildDfUniqueIndex, buildDfMultiIndex

# WARNING: UniProt is 1-based, open-end

class UniProtDataSet(SymEnum):
    SwissProt = 1
    TrEMBL = 2

class TransMatchStatus(SymEnum):
    no = 0
    maybe = 1
    yes = 2


def splitMetaList(val):
    """ split strings like: ENST00000369413.8|ENST00000528909.1"""
    if val == "":
        return []
    return val.split('|')

def splitDropVersion(idents):
    return [dropVersion(id) for id in splitMetaList(idents)]

def dropUniportIsoformModifier(acc):
    return acc.split('-')[0]

class UniProtMetaTbl:
    """reads swissprot.9606.tab, trembl.9606.tab"""
    def __init__(self, uniprotMetaTsv):
        self.df = pd.read_table(uniprotMetaTsv, keep_default_na=False)
        if "orgName" not in self.df.columns:
            raise Exception(f"column 'orgName' not found, is this a UnProt metadata TSV: {uniprotMetaTsv}")
        self.df.set_index('acc', inplace=True, drop=False, verify_integrity=True)
        self.geneNameIdx = buildDfMultiIndex(self.df, "geneName")

        # index by gene and transcript ids
        self.df['ensemblGeneIds'] = self.df.ensemblGene.apply(splitMetaList)
        self.df['ensemblGeneAccs'] = self.df.ensemblGene.apply(splitDropVersion)
        self.df['ensemblTransIds'] = self.df.ensemblTrans.apply(splitMetaList)
        self.df['ensemblTransAccs'] = self.df.ensemblTrans.apply(splitDropVersion)
        self.df['isoIds'] = self.df.isoIds.apply(splitMetaList)

        self.byGeneAccDf = self.df.explode('ensemblGeneAccs')
        self.byGeneAccDf.rename(columns={'ensemblGeneAccs': 'ensemblGeneAcc'}, inplace=True)
        self.byGeneAccDf.set_index('ensemblGeneAcc', inplace=True, drop=False, verify_integrity=False)

        self.byTranscriptAccDf = self.df.explode('ensemblTransAccs')
        self.byTranscriptAccDf.rename(columns={'ensemblTransAccs': 'ensemblTransAcc'}, inplace=True)
        self.byTranscriptAccDf.set_index('ensemblTransAcc', inplace=True, drop=False, verify_integrity=False)

        # with only one isoform: acc == mainIsoAcc, and isoIds are empty
        # for multiple isoforms: mainIsoAcc is canonical and isoIds are other isoforms.
        # Build an index by isoId which is
        self.df['allIsoIds'] = self.df.apply(lambda row: row.isoIds + [row.mainIsoAcc], axis=1)
        self.byIsoId = self.df.explode('allIsoIds')
        self.byIsoId.rename(columns={'allIsoIds': 'isoId'}, inplace=True)
        self.byIsoId.set_index('isoId', inplace=True, drop=False, verify_integrity=True)

    def getByAcc(self, acc):
        return self.df[self.df.index == acc].iloc[0]

    def getByIsoId(self, isoId):
        try:
            return self.byIsoId[self.byIsoId.index == isoId].iloc[0]
        except IndexError as ex:
            raise IndexError(f"isoId not found {isoId}") from ex

    def getTransMeta(self, transAcc):
        protMetas = self.byTranscriptAccDf[self.byTranscriptAccDf.index == transAcc]
        if len(protMetas) == 0:
            return None
        if len(protMetas) > 1:
            raise Exception(f"more than one protein for transcript {transAcc}, found: {protMetas.mainIsoAcc}")
        return protMetas.iloc[0]

    def getGeneAccMetas(self, geneAcc):
        protMetas = self.byGeneAccDf[self.byGeneAccDf.index == geneAcc]
        if len(protMetas) == 0:
            return None
        return protMetas

    def getGeneNameMetas(self, geneName):
        return self.geneNameIdx.get(geneName)

def uniprotParseAnnotId(annotId):
    "parse the annotation id into its parts: 'Q9BXI3#0' -> ('Q9BXI3' 0)"
    parts = annotId.rsplit('#', 1)
    if len(parts) != 2:
        raise Exception(f"invalid annotation id: '{annotId}'")
    return parts

def uniprotAnnotIdToAcc(annotId):
    "parse the annotation id into accession: 'Q9BXI3#0' -> 'Q9BXI3'"
    return uniprotParseAnnotId(annotId)[0]

class UniProtAnnotTbl:
    """reads swissprot.9606.annots.tab, trembl.9606.annots.tab"""
    def __init__(self, uniprotAnnotsTsv):
        self.df = pd.read_table(uniprotAnnotsTsv, keep_default_na=False)
        if "featType" not in self.df.columns:
            raise Exception("column 'featType' not found, is this a UnIprot annotations TSV: " + uniprotAnnotsTsv)

        # Add a unique, reproducible annotation id in the form mainIsoAcc#annotIdx, where the annotIdx is
        # the relative index of the annotation for that acc.
        self.df["annotId"] = self.df.mainIsoAcc.astype(str) + "#" + self.df.groupby("mainIsoAcc").cumcount().astype(str)

        self.annotIdIdx = buildDfUniqueIndex(self.df, "annotId")

    def getByAnnotId(self, annotId):
        annot = self.annotIdIdx.get(annotId)
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

def getAnnotColor(annot, dataSet):
    if dataSet == UniProtDataSet.TrEMBL:
        return TREMBLCOLOR
    color = commentColor.get(annot.comment, None)
    if color is None:
        color = featTypeColors.get(annot.featType, None)
    return color if color is not None else SWISSPCOLOR

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
        return ("unipMut", "Mutation")
    elif annot.featType in ["splice variant"]:
        return ("unipSplice", "Splice variant")
    elif annot.featType in ["strand", "helix", "coiled-coil region", "turn"]:
        return ("unipStruct", "Protein Structure")
    elif annot.featType in ["transmembrane region"]:
        # suggested by Regeneron: create four subtracks, for the subcell. localization indicators
        return ("unipLocTransMemb", "Transmembrane Domain")
    elif annot.comment in ["Extracellular"]:
        return ("unipLocExtra", "Extracellular Domain")
    elif annot.comment in ["Cytoplasmic"]:
        return ("unipLocCytopl", "Cytoplasmic Domains")
    elif annot.featType in ["signal peptide"]:
        return ("unipLocSignal", "Signal Peptide")
    elif annot.featType in ["chain"]:
        return ("unipChain", "Polypeptide Chains")
    elif annot.featType in ["repeat"]:
        return ("unipRepeat", "Repeat")
    elif annot.featType == "sequence conflict":
        return ("unipConflict", "Sequence Conflict")
    elif annot.featType in ["disulfide bond"]:
        return ("unipDisulfBond", "Disulfide Bond")
    elif annot.featType in ["domain", "zinc finger region", "topological domain"]:
        return ("unipDomain", "Domain")
    elif annot.featType in ["glycosylation site", "modified residue", "lipid moiety-binding region"]:
        return ("unipModif", "Amino Acid Modification")
    elif annot.featType in ["region of interest"]:
        return ("unipInterest", "Region of Interest")
    else:
        return ("unipOther", "Other Annotation")

def calcTransMatchStatus(uniprotMeta, transId):
    "determine match of transcript to GENCODE based on what transcripts are listed in metadata"
    if transId in uniprotMeta.ensemblTransIds:
        return TransMatchStatus.yes
    elif dropVersion(transId) in uniprotMeta.ensemblTransAccs:
        return TransMatchStatus.maybe
    else:
        return TransMatchStatus.no
