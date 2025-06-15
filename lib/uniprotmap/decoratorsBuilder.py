"""
Generate decorators using multiprocessing.  This calls an function to
generate decorators given PSLs in each subprocess.
"""
import multiprocessing as mp
import pipettor
from pycbio.sys import fileOps

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

def _worker(transAnnotMappings):
    """sub-process worker, if an error occurs an exception object is the returned.
    """
    if isinstance(_gAnnotationProcessor, Exception):
        return _gAnnotationProcessor
    try:
        return _gAnnotationProcessor.create(transAnnotMappings)
    except Exception as ex:
        ex2 = Exception("Worker failed")
        ex2.__cause__ = ex
        return ex

def _checkForWorkerFail(decoBeds):
    if isinstance(decoBeds, Exception):
        raise Exception("creation of decorator BEDs failed") from decoBeds

def _processSingle(annotationProcessorFactory,
                   transAnnotMappingReader, featTypeFunc,
                   decoBedFh, featTypes):
    # this is easier to debug without mp
    _workerInit(annotationProcessorFactory)
    for transAnnotMappings in transAnnotMappingReader:
        decoBeds = _worker(transAnnotMappings)
        _checkForWorkerFail(decoBeds)
        for decoBed in decoBeds:
            featTypes.add(featTypeFunc(decoBed))
            decoBed.write(decoBedFh)

def _processMulti(annotationProcessorFactory,
                  transAnnotMappingReader, featTypeFunc, nprocs,
                  decoBedFh, featTypes):
    with mp.Pool(processes=nprocs, initializer=_workerInit,
                 initargs=((annotationProcessorFactory,))) as pool:
        decoBedIters = pool.imap_unordered(_worker,
                                           transAnnotMappingReader)
        for decoBeds in decoBedIters:
            _checkForWorkerFail(decoBeds)
            for decoBed in decoBeds:
                featTypes.add(featTypeFunc(decoBed))
                decoBed.write(decoBedFh)

def _processMappings(annotationProcessorFactory, transAnnotMappingReader,
                     featTypeFunc, nprocs, decoBedFh, featTypes):
    # special case one process makes debugging & profiling easier
    if nprocs == 1:
        _processSingle(annotationProcessorFactory,
                       transAnnotMappingReader, featTypeFunc,
                       decoBedFh, featTypes)
    else:
        _processMulti(annotationProcessorFactory,
                      transAnnotMappingReader, featTypeFunc, nprocs,
                      decoBedFh, featTypes)


def buildDecorators(annotationProcessorFactory, transAnnotMappingReader,
                    featTypeFunc, annotDecoratorBedFile, nprocs):
    """
    Reads mapped annotation alignments and metadata for target transcripts, including
    those that did not align successfully. Yields TransAnnotMapping objects.

    Args:
        annot2GenomePslFile: Path to a PSL file mapping annotations to the genome.
        annot2GenomeRefTsv: Path to a TSV file with annotation metadata.
        annotLookupFunc: Callable that takes an annotation ID and returns its data record.
        transPslLookupFunc: Callable that takes a transcript ID and chromosome,
            and returns the corresponding alignment information.

    Yields:
        TransAnnotMapping: An object representing an annotation's genomic mapping.
        """
    featTypes = set()
    with fileOps.AtomicFileCreate(annotDecoratorBedFile) as tmpDecoBed:
        with pipettor.Popen(["sort"] + decoratorBedSortOpts, 'w', stdout=tmpDecoBed) as decoBedFh:
            _processMappings(annotationProcessorFactory, transAnnotMappingReader, featTypeFunc, nprocs,
                             decoBedFh, featTypes)
    return featTypes
