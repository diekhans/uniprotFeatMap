"""
Library of common functions and other definitions
"""

import sys
import os
from pycbio.sys import fileOps

def prMsg(msg):
    print(msg, file=sys.stderr, flush=True)

def dropVersion(ident):
    return ident.split('.')[0]


def mkPslMapCmd(inPsl, mapPsl, outPsl, *, swapMap=False, mapInfo=None,
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

class TmpOrSaveFile(str):
    """an file that might be saved for debugging purposes
    or discarded once used internally and then deleted"""
    def __new__(cls, interPrefix, suffix):
        if interPrefix is None:
            path = fileOps.tmpFileGet(suffix=suffix)
        else:
            path = interPrefix + suffix
        self = super(TmpOrSaveFile, cls).__new__(cls, path)
        self.isTmp = (interPrefix is None)
        return self

def cleanTmpFiles(*paths):
    """Remove each file is path is not None.  If it is a TmpOrSaveFile
    object, remove if it is not tmp"""
    for p in paths:
        if (p is not None) and (isinstance(p, str)) or (isinstance(p, TmpOrSaveFile) and p.isTmp):
            os.unlink(p)
