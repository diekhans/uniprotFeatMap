"""
common routing to support alignments on parasol
"""
from os import path as osp
import glob
import pipettor
from pycbio.distrib.parasol import Para
from pycbio.sys import fileOps
from protmap import conf, prMsg, getDoneFile


def queryGetSplitPrefix(queriesDir):
    return osp.join(queriesDir, "query")

def queryListSplitFas(queriesDir):
    return sorted(glob.glob(queryGetSplitPrefix(queriesDir) + "*"))

def querySplit(queryFa, queriesDir):
    prMsg("split query proteins")
    fileOps.ensureDir(queriesDir)
    pipettor.run(["faSplit", "about", queryFa, 2500, queryGetSplitPrefix(queriesDir)])

def makeJobFile(alignProg, queriesDir, targetDbFa, alignDir, alignBatchDir):
    jobFile = osp.join(alignBatchDir, "jobs.para")
    with open(jobFile, 'w') as fh:
        for queryFa in queryListSplitFas(queriesDir):
            outPsl = osp.join(alignDir, osp.basename(queryFa) + ".psl")
            print(alignProg, targetDbFa, queryFa, f"{{check out exists {outPsl}}}", file=fh)
    if osp.getsize(jobFile) == 0:
        raise Exception(f"empty job file create: {jobFile}")
    return jobFile

def runBatch(alignProg, queriesDir, targetDbFa, alignDir, alignBatchDir):
    prMsg("running alignment batch")
    fileOps.ensureDir(alignBatchDir)
    jobFile = makeJobFile(alignProg, queriesDir, targetDbFa, alignDir, alignBatchDir)
    para = Para(conf.paraHost, jobFile, paraDir=alignBatchDir)
    para.free()
    try:
        para.make()
    except pipettor.exceptions.ProcessException as ex:
        raise Exception(f"batch failed, correct problem, re-run with -batch={alignBatchDir}\n"
                        "then touch " + getDoneFile(alignDir)) from ex
