"""
Compares annotations on transcripts from different sources.
"""


class TransAnnot:
    "transcript annotations from UniProt, etc"
    def __init__(self, transId):
        self.transId = transId

class TargetTranscript:
    def __init__(self, transId):
        self.transId = transId
