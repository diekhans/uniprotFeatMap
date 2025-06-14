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
from uniprotmap import conf, prMsg
from uniprotmap.depends import runIfNotDone, runIfOutOfDate, getDoneFile

DEFAULT_QUERY_SPLIT_APPROX_SIZE = 25000

proteinTranscriptAlignJob = osp.normpath(osp.join(osp.dirname(__file__), "../../bin/proteinTranscriptAlignJob"))

class AlignError(Exception):
    pass

def updateCompoundFastaHeader(faRec):
    """convert >idA|idB|idC to >idA idB idC in a fasta Seq record"""
    if faRec.id.find('|') >= 0:
        idParts = faRec.id.split('|')
        faRec.id = idParts[0]
        desc = idParts[1:]
        if faRec.description != "":
            desc.append(faRec.description)
        faRec.description = ' '.join(desc)

##
# BLAST alignments
##
def _buildBlastTransIndex(transFa, workDir):
    logFile = osp.join(workDir, "formatdb.log")
    pipettor.run([osp.join(conf.blastDir, "formatdb"),
                  "-l", logFile, "-i", transFa, "-p", "F"])

##
# parasol batch alignments
##
def _makeJobFile(alignCmd, queriesDir, targetDbFa, alignDir, alignBatchDir):
    jobFile = osp.join(alignBatchDir, "jobs.para")
    with fileOps.opengz(jobFile, 'w') as fh:
        for queryFa in _queryListSplitFas(queriesDir):
            outPsl = osp.join(alignDir, osp.basename(queryFa) + ".psl")
            print(*alignCmd, targetDbFa, queryFa, f"{{check out exists {outPsl}}}", file=fh)
    if osp.getsize(jobFile) == 0:
        raise AlignError(f"empty job file create: {jobFile}")
    return jobFile

def _runBatch(alignCmdPre, queriesDir, targetDbFa, alignDir, alignBatchDir):
    "alignCmdPre is list of program and initial arguments"
    fileOps.ensureDir(alignBatchDir)
    jobFile = _makeJobFile(alignCmdPre, queriesDir, targetDbFa, alignDir, alignBatchDir)
    para = Para(paraHost=conf.paraHost, jobFile=jobFile, paraDir=alignBatchDir)
    para.clearSickNodes()
    para.freeBatch()
    try:
        para.make()
    except pipettor.exceptions.ProcessException as ex:
        raise AlignError(f"batch failed, correct problem, re-run with -batch={alignBatchDir}\n"
                         "then touch " + getDoneFile(alignDir)) from ex

##
# Alignment query setup
##
def _queryGetSplitPrefix(queriesDir):
    return osp.join(queriesDir, "query")

def _queryListSplitFas(queriesDir):
    return sorted(glob.glob(_queryGetSplitPrefix(queriesDir) + "*"))

def _queryFilterEditFasta(inFaFh, outFaFh, filterEditFunc):
    for faRec in SeqIO.parse(inFaFh, "fasta"):
        if (filterEditFunc is None) or filterEditFunc(faRec):
            SeqIO.write(faRec, outFaFh, "fasta")

def _queryBuildDb(queryFa, queriesDir, filterEditFunc, approxSize):
    """Split a query FASTA, If filterEditFunction is not none, it is passed the fasta record to
    check it should be included.  It can also update the FASTA record header if needed.
    """
    fileOps.ensureDir(queriesDir)
    # make sure there are no old files that could cause problems
    fileOps.rmFiles(*_queryListSplitFas(queriesDir))
    with fileOps.opengz(queryFa) as inFaFh:
        with pipettor.Popen(["faSplit", "about", "/dev/stdin", approxSize, _queryGetSplitPrefix(queriesDir)], 'w') as outFaFh:
            _queryFilterEditFasta(inFaFh, outFaFh, filterEditFunc)

##
# Alignment target setup
##
def _targetWriteMaskFastaRec(rec, outFh):
    """write a target transcript sequence where the CDS is upper case,
    hard-masking the """
    print(">" + rec.id, file=outFh)
    seq = re.sub('[a-z]', 'N', str(rec.seq))
    print(seq, file=outFh)

def _targetMakeUtrMaskedFasta(inFa, outFa, filterEditFunc):
    """
    Create FASTA with lower-case UTR hard-masked
    """
    with fileOps.opengz(inFa) as inFh:
        with fileOps.opengz(outFa, 'w') as outFh:
            for rec in SeqIO.parse(inFh, "fasta"):
                if (filterEditFunc is None) or filterEditFunc(rec):
                    _targetWriteMaskFastaRec(rec, outFh)

def _targetBuildDb(transFa, transDbFa, algo, filterEditFunc, targetDir):
    """build alignment target database for all transcripts were filterFunc(id) returns True"""
    fileOps.ensureDir(targetDir)
    fileOps.ensureFileDir(transDbFa)
    _targetMakeUtrMaskedFasta(transFa, transDbFa, filterEditFunc)
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

def _combinePairAligns(alignDir, prot2TransPslFile, filterFunc):
    "concatenate, filter, and sort by tName (transcript))"

    findSortCmd = (["find", alignDir, "-name", "*.fa.psl", "-print0"],
                   ["sort", "-k14,14", "-k10,10", "--files0-from=-"])
    with pipettor.Popen(findSortCmd, 'r') as inPslFh:
        with fileOps.opengz(prot2TransPslFile, 'w') as outPslFh:
            _processAlignedPsls(inPslFh, outPslFh, filterFunc)

##
# overall pipeline, parameterized with files and functions.
##
def proteinTranscriptAlign(protFa, transFa, prot2TransPslFile, algo, workDir, *,
                           queryFaEditFilterFunc=None, targetFaEditFilterFunc=None, alignFilterFunc=None,
                           querySplitSize=DEFAULT_QUERY_SPLIT_APPROX_SIZE):
    targetDir = osp.join(workDir, "transDb")
    transDbFa = osp.join(targetDir, "transDb.fa")
    with runIfNotDone(targetDir, depends=transFa) as do:
        if do:
            prMsg("building target transcript database")
            _targetBuildDb(transFa, transDbFa, algo, targetFaEditFilterFunc, targetDir)

    queryDir = osp.join(workDir, "queryDir")
    with runIfNotDone(queryDir, depends=protFa) as do:
        if do:
            prMsg("split proteins")
            _queryBuildDb(protFa, queryDir, queryFaEditFilterFunc, querySplitSize)

    alignDir = osp.join(workDir, "aligns")
    with runIfNotDone(alignDir, doneDepends=[targetDir, queryDir]) as do:
        if do:
            prMsg("running alignment batch")
            _runBatch([proteinTranscriptAlignJob, algo], queryDir, transDbFa, alignDir,
                      osp.join(workDir, "batch"))

    with runIfOutOfDate(prot2TransPslFile, doneDepends=alignDir) as do:
        if do:
            prMsg("combining alignments")
            with fileOps.AtomicFileCreate(prot2TransPslFile) as tmpPslFile:
                _combinePairAligns(alignDir, tmpPslFile, alignFilterFunc)
    prMsg("finished")
