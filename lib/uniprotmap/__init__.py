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
