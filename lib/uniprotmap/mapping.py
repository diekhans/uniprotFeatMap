##
# pslMap related functions
##
from collections import defaultdict
from pycbio.tsv import TsvReader, strOrNoneType, intOrNoneType

def pslMapMkCmd(inPsl, mapPsl, outPsl, *, swapMap=False, mapInfo=None,
                interPrefix=None, interMid=None, chainMapFile=False):
    """Build pslMap command, possibly saving intermediate files for debugging.
    Intermediate files are save if interPrefix is not None, and include
    interMid.  interMid is ignored if interPrefix is None.
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

    if inter is not None:
        if mapInfo is None:
            cmd1.append(f"-mapInfo={inter}.mapInfo.tsv")
        cmd2 = ["tee", f"{inter}.psl"]
        cmds.append(cmd2)
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
        self.bySrcQName = defaultdict(list)
        self.byMappedTName = defaultdict(list)
        for row in TsvReader(mapInfoTsv, typeMap=_mapInfoTypeMap):
            self.append(row)
            self.bySrcQName[row.srcQName].append(row)
            if row.mappedTName is not None:
                self.byMappedTName[row.mappedTName].append(row)
        self.bySrcQName.default_factory = None
        self.byMappedTName.default_factory = None
