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

def annotIdFmt(protAcc, annotIdx):
    return f"{protAcc}|{annotIdx}"

def annotIdParse(annotId):
    parts = annotId.split('|')
    try:
        if len(parts) != 2:
            raise ValueError(f"expected 'protAcc|annotIdx', got '{annotId}'")
        return parts[0], int(parts[1])
    except Exception as ex:   # catch bad int
        raise ValueError(f"invalid annotation id: '{annotId}'") from ex

def annotIdToProtAcc(annotId):
    return annotIdParse(annotId)[0]

def annotMapIdFmt(annotId, mapIdx):
    return f"{annotId}|{mapIdx}"

def annotMapIdParse(annotMapId):
    parts = annotMapId.split('|')
    try:
        if len(parts) != 3:
            raise ValueError(f"expected 'protAcc|annotIdx|mapIdx', got '{annotMapId}'")
        return parts[0], int(parts[1]), int(parts[2])
    except Exception as ex:   # catch bad int
        raise ValueError(f"invalid annotation map id: '{annotMapId}'") from ex

def annotMapIdToAnnotId(annotMapId):
    protAcc, annotIdx, _ = annotMapIdParse(annotMapId)
    return annotIdFmt(protAcc, annotIdx)

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
        if (p is not None) and p.isTmp:
            os.unlink(p)
