"""
Support for creating UniProt decorators
"""

from pycbio.sys.svgcolors import SvgColors
from pycbio.sys.color import Color
from pycbio.sys.symEnum import SymEnum, auto
from uniprotmap import dropVersion
from uniprotmap.uniprot import UniProtCategory, UniProtDataSet, TransCategory
from uniprotmap.mappingAnalysis import FeatureIndelType
from pycbio.hgdata.bed import encodeRow
from pycbio.hgdata.decoration import Decoration

def makeColorDesc(color, info):
    colorName = SvgColors.getClosestName(color)
    colorRgb = color.toRgb8Str() if color.alpha is None else color.toRgba8Str()
    return f"{info:<36} {colorName:<14} {colorRgb}"

######
# The below code is derived from kent/src/hg/utils/otto/uniprot/doUniprot; if
# this converges, we can turn this into a library.
######

def isMutagenesis(annot):
    return annot.featType == "mutagenesis site"

def isVariant(annot):
    return annot.featType == "sequence variant"


###
# Color logic.  This mostly matches UniProt genomic track, however colors for
# overlay are not visible with track colors, so these are adjusted.
###

def _mkcolor(r, g, b, a=None):
    return Color.fromRgb8(r, g, b, a)

UNIPROT_CANON_ISO_OUTLINE_COLOR = _mkcolor(0, 191, 255)    # deepskyblue
UNIPROT_NONCANON_ISO_OUTLINE_COLOR = _mkcolor(0, 0, 0)     # black
FEAT_INSERTION_COLOR = _mkcolor(255, 0, 0)                 # red
FEAT_DELETION_COLOR = _mkcolor(255, 0, 0)                  # red

SWISSPCOLOR = _mkcolor(34, 139, 34)    # forestgreen, genomic tracks are 12,12,120 dark blue
TREMBLCOLOR = _mkcolor(143, 188, 143)  # darkseagreen, genomic tracks are 0,150,250 light blue

# mapping of annotations columns to colors
featTypeColors = {
    "modified residue": _mkcolor(255, 215, 0),            # gold
    "glycosylation site": _mkcolor(0, 100, 100),          # teal
    "disulfide bond": _mkcolor(100, 100, 100),            # dimgray
    "topological domain": _mkcolor(100, 0, 0),            # maroon
    "zinc finger region": _mkcolor(100, 100, 0),          # olive
    "transmembrane region": _mkcolor(0, 150, 0),          # green
    "signal peptide": _mkcolor(255, 0, 150),              # deeppink
    "lipid moiety-binding region": _mkcolor(12, 12, 120)  # navy
}
shortFeatTypeColors = {
    "phos": _mkcolor(200, 200, 0),        # goldenrod
}
commentColor = {
    "Extracellular": _mkcolor(0, 110, 180),  # darkcyan
    "Cytoplasmic": _mkcolor(255, 150, 0),    # darkorange
}

# TrEMBL categories to provide special coloring
_tremblColorCategories = frozenset((
    UniProtCategory.LocSignal,
    UniProtCategory.LocExtra,
    UniProtCategory.LocTransMemb,
    UniProtCategory.LocCytopl,
))

def getFeatTypeColor(annot):
    "by feature or short feature type"
    color = featTypeColors.get(annot.featType, None)
    if color is None:
        color = shortFeatTypeColors.get(annot.shortFeatType, None)
    return color

def getAnnotColorSwissport(annot):
    color = commentColor.get(annot.comment, None)
    if color is None:
        color = getFeatTypeColor(annot)
    return color if color is not None else SWISSPCOLOR

def getAnnotColorTrembl(annot):
    if getAnnotCategory(annot)[0] in _tremblColorCategories:
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
    if annot.shortFeatType == "phos":
        return "Phosphorylation"
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
    elif annot.featType in ["glycosylation site", "phosphorylation site", "modified residue", "lipid moiety-binding region"]:
        return (UniProtCategory.Modif, "Amino Acid Modification")
    elif annot.featType in ["region of interest"]:
        return (UniProtCategory.Interest, "Region of Interest")
    else:
        return (UniProtCategory.Other, "Other Annotation")

def calcTransCategory(uniprotMeta, transId):
    "determine match of transcript to based on what transcripts are listed in metadata"
    if transId in uniprotMeta.ensemblTransIds:
        return TransCategory.canonical
    if dropVersion(transId) in uniprotMeta.ensemblTransAccs:
        return TransCategory.canonical_diff_version
    else:
        return TransCategory.noncanonical

def getColorUses():
    "list of tuples of (color, use_description)"
    uses = [(SWISSPCOLOR, "SwissProt base"),
            (TREMBLCOLOR, "TrEMBL base")]
    for info, color in featTypeColors.items():
        uses.append((color, info))
    for info, color in shortFeatTypeColors.items():
        uses.append((color, info))
    for info, color in commentColor.items():
        uses.append((color, info))
    return uses

# sort command options for sorting decorator BED. Include fields for test
# reproducible and locality

class AnnotType(SymEnum):
    feature = auto()
    disruption = auto()

class FeatStatus(SymEnum):
    complete = auto()
    disrupted = auto()
    insertion = auto()
    deletion = auto()

class UniprotDecoration(Decoration):
    """Holds a decoration BED + extra data.  This matches etc/uniprotDecoration.as. """

    # Keep in schema order:
    __slots__ = ("annotType", "dataSet", "uniprotAcc", "transCategory", "canonTransId",
                 "featStatus", "category", "categoryName", "description", "shortFeatType", "featType",
                 "shortName", "longName", "comment", "disease")

    def __init__(self, chrom, bedBlocks, name, strand, color,
                 itemName, itemStart, itemEnd, glyph, fillColor, *, annotType,
                 dataSet, uniprotAcc, transCategory, canonTransId, featStatus,
                 category, categoryName, description, shortFeatType, featType,
                 shortName, longName, comment, disease):
        # sanity check
        assert isinstance(dataSet, UniProtDataSet)
        assert isinstance(transCategory, TransCategory)
        assert isinstance(annotType, AnnotType)
        assert isinstance(featStatus, (FeatStatus, FeatureIndelType))
        assert isinstance(category, UniProtCategory)

        super(UniprotDecoration, self).__init__(chrom, bedBlocks[0].start, bedBlocks[-1].end, name,
                                                itemName, itemStart, itemEnd,
                                                strand=strand, itemRgb=color, blocks=bedBlocks,
                                                glyph=glyph, fillColor=fillColor)

        self.annotType = annotType
        self.dataSet = dataSet
        self.uniprotAcc = uniprotAcc
        self.transCategory = transCategory
        self.canonTransId = canonTransId
        self.featStatus = featStatus
        self.category = category
        self.categoryName = categoryName
        self.description = description
        self.shortFeatType = shortFeatType
        self.featType = featType
        self.shortName = shortName
        self.longName = longName
        self.comment = comment
        self.disease = disease

    def toRow(self):
        return super().toRow() + encodeRow([getattr(self, f) for f in self.__slots__])
