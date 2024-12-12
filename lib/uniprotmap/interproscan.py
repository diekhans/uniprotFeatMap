"""
Access to interproscan results in JSON format.
"""
import os.path as osp
import json
import pickle
from pycbio.sys import fileOps
from pycbio.sys.objDict import ObjDict

class InterproError(Exception):
    pass

class InterproResults:
    """Interproscan analysis results from JSON output"""
    def __init__(self):
        self.byProtId = {}

    def add(self, result):
        self.byProtId[result.xref.id] = result

def _interproPklFile(interproJson, cacheDir):
    """removes compression extension and .json, if present, and return .pkl.gz
    in cache dir"""
    inBase = osp.basename(fileOps.compressBaseName(interproJson))
    return osp.join(cacheDir, inBase + ".pkl.gz")

def _interproReadJson(interproJson):
    interproResults = InterproResults()
    with fileOps.opengz(interproJson) as fh:
        for result in json.load(fh, object_pairs_hook=ObjDict):
            interproResults.add(result)
    return interproResults

def _interproPklCurrent(interproJson, interproPkl):
    return osp.exists(interproPkl) and (osp.getmtime(interproJson) <= osp.getmtime(interproPkl))

def _interproSavePkl(interproResults, interproPkl):
    with fileOps.opengz(interproPkl, 'wb') as fh:
        pickle.dump(interproResults, fh)

def _interproLoadPkl(interproPkl):
    with fileOps.opengz(interproPkl, 'rb') as fh:
        return pickle.load(fh)


def interproResultsLoad(interproJson, *, cacheDir=None):
    """Load interproscan results. If cache_dir is not None, it is used to
    pickle the InterproscanResults object, with replacing the .json[.gz], with .pkl.gz.
    If pickle file does not exist or is out-of-date, it is creates """

    if cacheDir is None:
        return _interproReadJson(interproJson)
    else:
        interproPkl = _interproPklFile(interproJson, cacheDir)
        if _interproPklCurrent(interproJson, interproPkl):
            return _interproLoadPkl(interproPkl)
        else:
            interproResults = _interproReadJson(interproJson)
            _interproSavePkl(interproResults, interproPkl)
            return interproResults
