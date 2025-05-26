"""
Compares annotations on transcripts from different sources.
"""


class AnnotSet:
    def __init__(self, setName):
        set.setName = setName

class AnnotMap:
    """Mappings of an annotation to target genome"""
    __slots__ = ('annotId', 'mapRanges')

    def __init__(self, annotId, mapRanges):
        self.annotId = annotId
        self.mapRanges = mapRanges


class TransAnnots:
    "transcript annotations from UniProt, etc"
    def __init__(self, annotSet, transId):
        self.annotSet = annotSet
        self.transId = transId

class TargetTranscript:
    def __init__(self, transId):
        self.transId = transId


# def _mapTransAnnots(transGenomePsl, annotTransPsls):
#    for annotTransPsl in annotTransPsls:
