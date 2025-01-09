"""
Generate decorators using multiprocessing.  This calls an function to
generate decorators given PSLs in each subprocess.
"""
import multiprocessing as mp
import pipettor
from pycbio.sys import fileOps
from pycbio.hgdata.psl import PslReader

# chrom,  chromStart, chromEnd, decoratedItem, name, dataset
decoratorBedSortOpts = ["-k1,1", "-k2,2n", "-k3,3n", "-k13,13", "-k4,4n", "-k17,17"]

##
# This holds the forked process-global instance of AnnotationProcessor
# If an error occurs during initialization, it is set to the Exception
##
_gAnnotationProcessor = None

def _workerInit(annotationProcessorFactory):
    "sub-process setup, will set global annotation process or store an exception."
    global _gAnnotationProcessor
    try:
        _gAnnotationProcessor = annotationProcessorFactory()
        assert _gAnnotationProcessor is not None
    except Exception as ex:
        _gAnnotationProcessor = Exception("Pool initialization failed")
        _gAnnotationProcessor.__cause__ = ex

def _worker(alignBatch):
    """sub-process worker, if an error occurs an exception object is the returned.
    """
    if isinstance(_gAnnotationProcessor, Exception):
        return gAnnotationProcessor
    try:
        decoBeds = []
        for alignIdx, annotPsl in alignBatch:
            beds = _gAnnotationProcessor.create(alignIdx, annotPsl)
            if beds is not None:
                decoBeds.extend(beds)
        return decoBeds
    except Exception as ex:
        ex2 = Exception("Worker failed")
        ex2.__cause__ = ex
        return ex

def _annotGenomeBatchReader(annotGenomePslFile, batchSize):
    """return tuples of (alignIdx, annotPsl)"""
    alignIdx = 0
    alignBatch = []

    for annotPsl in PslReader(annotGenomePslFile):
        alignBatch.append((alignIdx, annotPsl))
        alignIdx += 1
        if len(alignBatch) >= batchSize:
            yield alignBatch
            alignBatch = []
    if len(alignBatch) > 0:
        yield alignBatch

def _checkForWorkerFail(decoBeds):
    if isinstance(decoBeds, Exception):
        raise Exception("creation of decorator BEDs failed") from decoBeds

def _writeDecoratorBeds(featTypeFunc, decoBedIters, annotDecoratorBedFile):
    """write results from multiprocessing task as they come back via an iter
    over lists of beds, or maybe a list of lists. This also checks for
    exceptions in the list, since these are combine back asynchronously, we
    have to way until actually writing."""
    featTypes = set()
    with fileOps.AtomicFileCreate(annotDecoratorBedFile) as tmpDecoBed:
        with pipettor.Popen(["sort"] + decoratorBedSortOpts, 'w', stdout=tmpDecoBed) as decoBedFh:
            for decoBeds in decoBedIters:
                _checkForWorkerFail(decoBeds)
                for decoBed in decoBeds:
                    decoBed.write(decoBedFh)
                    featTypes.add(featTypeFunc(decoBed))
    return featTypes

def processSingle(annotationProcessorFactory, featTypeFunc,
                  annotGenomePslFile, annotTransRefTsv, transGenomePsl,
                  annotDecoratorBedFile, batchSize):
    # this is easier to debug without mp
    _workerInit(annotationProcessorFactory)
    decoBedsList = []
    for alignBatch in _annotGenomeBatchReader(annotGenomePslFile, batchSize):
        decoBeds = _worker(alignBatch)
        _checkForWorkerFail(decoBeds)
        decoBedsList.append(decoBeds)

    featTypes = _writeDecoratorBeds(featTypeFunc, decoBedsList, annotDecoratorBedFile)
    return featTypes

def processMulti(annotationProcessorFactory, featTypeFunc,
                 annotGenomePslFile, annotTransRefTsv, transGenomePsl,
                 annotDecoratorBedFile, nprocs, batchSize):
    with mp.Pool(processes=nprocs, initializer=_workerInit,
                 initargs=((annotationProcessorFactory,))) as pool:
        decoBedIters = pool.imap_unordered(_worker,
                                           _annotGenomeBatchReader(annotGenomePslFile, batchSize))

        featTypes = _writeDecoratorBeds(featTypeFunc, decoBedIters, annotDecoratorBedFile)
    return featTypes

def buildDecorators(annotationProcessorFactory, featTypeFunc,
                    annotGenomePslFile, annotTransRefTsv, transGenomePsl,
                    annotDecoratorBedFile, nprocs, batchSize):
    """
    Used to build decorators using multiple processes.  Special handling
    for nprocs=1 to build in current process to make debugging easier.

    - annotationProcessorFactory function, usually a partial with arguments,
      creates an object convert PSLs to decorators.  It should have one function
        annotationProcessor.create(alignIdx, annotPsl)
      that either returns a *list* of decorator BED records or None.
    - featTypeFunc - a function that given a decorate BED returns an tuple
      that describes features used in the decorates to use in building filters.
    """
    # special case one process makes profiling easier
    if nprocs == 1:
        featTypes = processSingle(annotationProcessorFactory, featTypeFunc,
                                  annotGenomePslFile, annotTransRefTsv, transGenomePsl,
                                  annotDecoratorBedFile, batchSize)
    else:
        featTypes = processMulti(annotationProcessorFactory, featTypeFunc,
                                 annotGenomePslFile, annotTransRefTsv, transGenomePsl,
                                 annotDecoratorBedFile, nprocs, batchSize)
    return featTypes
