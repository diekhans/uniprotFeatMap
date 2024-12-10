"""
common routing to support alignments on parasol
"""
from os import path as osp
import re
import glob
import pipettor
from pycbio.distrib.parasol import Para
from pycbio.sys import fileOps
from pycbio.hgdata.psl import PslReader
from Bio import SeqIO
from uniprotmap import conf
from uniprotmap.depends import getDoneFile

DEFAULT_QUERY_SPLIT_APPROX_SIZE = 25000

class AlignError(Exception):
    pass

##
# BLAST alignments
##
def _buildBlastTransIndex(transFa, workDir):
    logFile = osp.join(workDir, "formatdb.log")
    pipettor.run([osp.join(conf.blastDir, "formatdb"),
                  "-l", logFile, "-i", transFa, "-p", "F"])

##
# Alignment query setup
##
def queryGetSplitPrefix(queriesDir):
    return osp.join(queriesDir, "query")

def queryListSplitFas(queriesDir):
    return sorted(glob.glob(queryGetSplitPrefix(queriesDir) + "*"))

def _querySplitWriter(inFaFh, outFaFh, filterEditFunc):
    for faRec in SeqIO.parse(inFaFh, "fasta"):
        if (filterEditFunc is None) or filterEditFunc(faRec):
            SeqIO.write(faRec, outFaFh, "fasta")

def queryBuildDb(queryFa, queriesDir, *, filterEditFunc=None, approxSize=DEFAULT_QUERY_SPLIT_APPROX_SIZE):
    """Split a query FASTA, If filterEditFunction is not none, it is passed the fasta record to
    check it should be included.  It can also update the FASTA record header if needed.
    """
    fileOps.ensureDir(queriesDir)
    # make sure there are no old files that could cause problems
    fileOps.rmFiles(*queryListSplitFas(queriesDir))
    with fileOps.opengz(queryFa) as inFaFh:
        with pipettor.Popen(["faSplit", "about", "/dev/stdin", approxSize, queryGetSplitPrefix(queriesDir)], 'w') as outFaFh:
            _querySplitWriter(inFaFh, outFaFh, filterEditFunc)

##
# Alignment target setup
##
def _targetWriteMaskFastaRec(rec, outFh):
    """write a target transcript sequence where the CDS is upper case,
    hard-masking the """
    print(">" + rec.id, file=outFh)
    seq = re.sub('[a-z]', 'N', str(rec.seq))
    print(seq, file=outFh)

def _targetMakeUtrMaskedFasta(filterFunc, inFa, outFa):
    """
    Create FASTA with lower-case UTR hard-masked
    """
    with fileOps.opengz(inFa) as inFh:
        with fileOps.opengz(outFa, 'w') as outFh:
            for rec in SeqIO.parse(inFh, "fasta"):
                if filterFunc(rec.id):
                    _targetWriteMaskFastaRec(rec, outFh)

def targetBuildDb(filterFunc, transFa, transDbFa, algo, targetDir):
    """build alignment target database for all transcripts were filterFunc(id) returns True"""
    fileOps.ensureDir(targetDir)
    fileOps.ensureFileDir(transDbFa)
    _targetMakeUtrMaskedFasta(filterFunc, transFa, transDbFa)
    if algo == "blast":
        _buildBlastTransIndex(transDbFa, targetDir)

##
# alignment results collection
##
def _queryTargetPairPslReader(inPslFh):
    """read batches with same set of name query and target names Input mushed
    be sorted by (target, query)."""
    # Batch program has already discards PSLs that are not ++ alignments,
    pairedPsls = []
    for psl in PslReader(inPslFh):
        if len(pairedPsls) == 0:
            pairedPsls.append(psl)
        elif ((psl.qName == pairedPsls[0].qName) and
              (psl.tName == pairedPsls[0].tName)):
            pairedPsls.append(psl)
        else:
            yield pairedPsls
            pairedPsls = [psl]
    if len(pairedPsls) > 0:
        yield pairedPsls

def _selectPairedPsls(pairedPsls, filterFunc):
    # all have same query/target
    if not filterFunc(pairedPsls[0]):
        return None
    if len(pairedPsls) > 1:
        # most query aligned wins
        pairedPsls.sort(key=lambda p: p.queryAligned(), reverse=True)
    return pairedPsls[0]

def _processAlignedPsls(inPslFh, outPslFh, filterFunc):
    for pairedPsls in _queryTargetPairPslReader(inPslFh):
        psl = _selectPairedPsls(pairedPsls, filterFunc)
        if psl is not None:
            psl.write(outPslFh)

def combinePairAligns(alignDir, protTransPsl, filterFunc):
    "concatenate, filter, and sort by tName (transcript))"

    findSortCmd = (["find", alignDir, "-name", "*.fa.psl", "-print0"],
                   ["sort", "-k14,14", "-k10,10", "--files0-from=-"])
    with pipettor.Popen(findSortCmd, 'r') as inPslFh:
        with fileOps.opengz(protTransPsl, 'w') as outPslFh:
            _processAlignedPsls(inPslFh, outPslFh, filterFunc)

##
# parasol batch alignments
##
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
    para = Para(conf.paraHost, jobFile=jobFile, paraDir=alignBatchDir, retries=2)
    para.free()
    try:
        para.make()
    except pipettor.exceptions.ProcessException as ex:
        raise AlignError(f"batch failed, correct problem, re-run with -batch={alignBatchDir}\n"
                         "then touch " + getDoneFile(alignDir)) from ex
