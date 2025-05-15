##
# Common function for mapping of annotations to the genome using pslMap
# and related
##
from collections import defaultdict
import pipettor
from pycbio.tsv import TsvReader, strOrNoneType, intOrNoneType
from pycbio.hgdata.psl import Psl, PslBlock, PslReader

def pslMapMkCmd(inPsl, mapPsl, outPsl, *, swapMap=False, outPslCopy=None, mapInfo=None,
                interPrefix=None, interMid=None, chainMapFile=False):
    """Build pslMap command, possibly saving intermediate files for debugging.
    Intermediate files are save if interPrefix is not None, and include
    interMid.  interMid is ignored if interPrefix is None.
    The outPslCopy is used to save output when main output is going to
    stdout.
    """
    assert (interPrefix is None) or (interMid is not None)
    inter = None if interPrefix is None else f"{interPrefix}{interMid}"

    cmd1 = ["pslMap", "-tsv", "-check", "-inType=na_na", "-mapType=na_na"]
    cmds = [cmd1]
    if swapMap:
        cmd1.append("-swapMap")
    if chainMapFile:
        cmd1.append("-chainMapFile")
    if mapInfo is not None:
        cmd1.append(f"-mapInfo={mapInfo}")
    cmd1.extend((inPsl, mapPsl, outPsl))

    if outPslCopy is not None:
        cmds.append(["tee", outPslCopy])
    if inter is not None:
        if mapInfo is None:
            cmd1.append(f"-mapInfo={inter}.mapInfo.tsv")
        cmds.append(["tee", f"{inter}.psl"])
    return cmds

_mapInfoTypeMap = {
    'srcQStart': int,
    'srcQEnd': int,
    'srcQSize': int,
    'srcTStart': int,
    'srcTEnd': int,
    'srcAligned': int,
    'mappingQName': strOrNoneType,
    'mappingQStart': intOrNoneType,
    'mappingQEnd': intOrNoneType,
    'mappingTName': strOrNoneType,
    'mappingTStart': intOrNoneType,
    'mappingTEnd': intOrNoneType,
    'mappingId': intOrNoneType,
    'mappedQName': strOrNoneType,
    'mappedQStart': intOrNoneType,
    'mappedQEnd': intOrNoneType,
    'mappedTName': strOrNoneType,
    'mappedTStart': intOrNoneType,
    'mappedTEnd': intOrNoneType,
    'mappedAligned': intOrNoneType,
    'qStartTrunc': intOrNoneType,
    'qEndTrunc': intOrNoneType,
    'mappedPslLine': intOrNoneType,
}

class PslMapInfoTbl(list):
    """read and index pslMap -mapInfo files"""
    def __init__(self, mapInfoTsv):
        self.bySrcTName = defaultdict(list)
        # use (mappingQName, mappedTName) to handle PAR
        self.byMappingQMappedTNames = defaultdict(list)
        self.byMappedTName = defaultdict(list)
        for row in TsvReader(mapInfoTsv, typeMap=_mapInfoTypeMap):
            self._load_row(row)
        self.bySrcTName.default_factory = None
        self.byMappingQMappedTNames.default_factory = None
        self.byMappedTName.default_factory = None

    def _load_row(self, row):
        self.append(row)
        self.bySrcTName[row.srcTName].append(row)
        if row.mappedTName is not None:
            self.byMappingQMappedTNames[(row.mappingQName, row.mappedTName)].append(row)
            self.byMappedTName[row.mappedTName].append(row)

def getQuerySizes(pslFile):
    """Get the sizes of all of queries in a PSL file."""
    querySizes = {}
    for psl in PslReader(pslFile):
        querySizes[psl.qName] = psl.qSize
    return querySizes

def createAnnotToProteinCdsPsl(annotId, annotStartOff, annotEndOff,
                               protId, protCdsSize):
    """Create a PSL alignment of a protein annotation to it's protein.
    The size and coordinates should be converted to bases (3x AA).
    Results are in nucleic acid (CDS) coordinates."""
    annotSize = annotEndOff - annotStartOff
    psl = Psl(qName=annotId, qSize=annotSize, qStart=0, qEnd=annotSize,
              tName=protId, tSize=protCdsSize, tStart=annotStartOff, tEnd=annotEndOff,
              strand='+')
    psl.addBlock(PslBlock(qStart=0, tStart=annotStartOff, size=annotSize))
    psl.computeCounts()
    return psl

def pslMapAnnots(annotCanonPsl, prot2TransPairedPsl, trans2GenomePsl,
                 annotGenomeMapInfoTsv, annot2GenomePslFh, *, annot2TransPsl=None,
                 interPrefix=None, xspeciesTrans2TransPsl=None, xspeciesTransMapInfoTsv=None):
    # all of these pslMaps are NA/NA -> NA/NA -> NA/NA
    cmds = []

    # annotation on canonical transcripts to all transcripts alignment to canonical
    cmds += pslMapMkCmd(annotCanonPsl, prot2TransPairedPsl, "/dev/stdout",
                        outPslCopy=annot2TransPsl, interPrefix=interPrefix, interMid="annotTrans")
    if xspeciesTrans2TransPsl is not None:
        # transcript to other species transcript mapping
        cmds += pslMapMkCmd("/dev/stdin", xspeciesTrans2TransPsl, "/dev/stdout",
                            mapInfo=xspeciesTransMapInfoTsv,
                            interPrefix=interPrefix, interMid="xspeciesAnnotTrans")

    # per-transcript annotation to genome mapping
    cmds += pslMapMkCmd("/dev/stdin", trans2GenomePsl, "/dev/stdout", mapInfo=annotGenomeMapInfoTsv)

    pipettor.run(cmds, stdout=annot2GenomePslFh)
