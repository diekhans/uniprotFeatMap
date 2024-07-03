"""
common routing to support alignments on parasol
"""
from os import path as osp
import glob
import pipettor
from pycbio.distrib.parasol import Para
from pycbio.sys import fileOps
from uniprotmap import conf
from uniprotmap.depends import getDoneFile

class AlignError(Exception):
    pass

def buildBlastTransIndex(transFa, workDir):
    logFile = osp.join(workDir, "formatdb.log")
    pipettor.run([osp.join(conf.blastDir, "formatdb"),
                  "-l", logFile, "-i", transFa, "-p", "F"])

def queryGetSplitPrefix(queriesDir):
    return osp.join(queriesDir, "query")

def queryListSplitFas(queriesDir):
    return sorted(glob.glob(queryGetSplitPrefix(queriesDir) + "*"))

def _querySplitWriter(inFaFh, outFaFh, filterFunc):
    incl = False
    for line in inFaFh:
        if line.startswith('>'):
            incl = filterFunc(line)
        if incl:
            outFaFh.write(line)

def querySplit(queryFa, queriesDir, *, filterFunc, approxSize=25000):
    fileOps.ensureDir(queriesDir)
    # make sure there are no old files that could corrupt
    fileOps.rmFiles(*queryListSplitFas(queriesDir))
    with fileOps.opengz(queryFa) as inFaFh:
        with pipettor.Popen(["faSplit", "about", "/dev/stdin", approxSize, queryGetSplitPrefix(queriesDir)], 'w') as outFaFh:
            _querySplitWriter(inFaFh, outFaFh, filterFunc)

def makeJobFile(alignCmd, queriesDir, targetDbFa, alignDir, alignBatchDir):
    jobFile = osp.join(alignBatchDir, "jobs.para")
    with fileOps.opengz(jobFile, 'w') as fh:
        for queryFa in queryListSplitFas(queriesDir):
            outPsl = osp.join(alignDir, osp.basename(queryFa) + ".psl")
            print(*alignCmd, targetDbFa, queryFa, f"{{check out exists {outPsl}}}", file=fh)
    if osp.getsize(jobFile) == 0:
        raise AlignError(f"empty job file create: {jobFile}")
    return jobFile

def runBatch(alignCmdPre, queriesDir, targetDbFa, alignDir, alignBatchDir):
    "alignCmdPre is list of program and initial arguments"
    fileOps.ensureDir(alignBatchDir)
    jobFile = makeJobFile(alignCmdPre, queriesDir, targetDbFa, alignDir, alignBatchDir)
    para = Para(conf.paraHost, jobFile=jobFile, paraDir=alignBatchDir, retries=0)
    para.free()
    try:
        para.make()
    except pipettor.exceptions.ProcessException as ex:
        raise AlignError(f"batch failed, correct problem, re-run with -batch={alignBatchDir}\n"
                         "then touch " + getDoneFile(alignDir)) from ex
