"""
common routing to support alignments on parasol
"""
from os import path as osp
import glob
import pipettor
from pycbio.distrib.parasol import Para
from pycbio.sys import fileOps
from protmap import conf, prMsg, getDoneFile


def queryGetSplitPrefix(querySplitDir):
    return osp.join(querySplitDir, "query")

def queryListSplitFas(querySplitDir):
    return glob.glob(queryGetSplitPrefix(querySplitDir) + "*")

def querySplit(queryFa, querySplitDir):
    prMsg("split query proteins")
    fileOps.ensureDir(querySplitDir)
    pipettor.run(["faSplit", "about", queryFa, 2500, queryGetSplitPrefix(querySplitDir)])

def makeJobFile(bindir, querySplitDir, transDbFa, alignDir, workDir):
    jobFile = osp.join(workDir, "jobs.para")
    jobCmd = osp.join(bindir, "proteinTranscriptAlignJob")
    with open(jobFile, 'w') as fh:
        for queryFa in queryListSplitFas(querySplitDir):
            outPsl = osp.join(alignDir, osp.basename(queryFa) + ".psl")
            print(jobCmd, transDbFa, queryFa, f"{{check out exists {outPsl}}}", file=fh)
    if osp.getsize(jobFile) == 0:
        raise Exception(f"empty job file create: {jobFile}")
    return jobFile

def runBatch(bindir, querySplitDir, transDbFa, alignDir, workDir):
    prMsg("running alignment batch")
    jobFile = makeJobFile(bindir, querySplitDir, transDbFa, alignDir, workDir)
    paraDir = osp.abspath(osp.join(workDir, "batch"))
    para = Para(conf.paraHost, jobFile, paraDir=paraDir)
    para.free()
    try:
        para.make()
    except pipettor.exceptions.ProcessException as ex:
        raise Exception(f"batch failed, correct problem, re-run with -batch={paraDir}\n"
                        "then touch " + getDoneFile(alignDir)) from ex
