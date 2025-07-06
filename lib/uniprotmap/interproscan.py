"""
Access to interproscan results in JSON format.
"""
from collections import defaultdict
from pycbio.tsv import TsvReader, TsvRow
from uniprotmap import annotIdFmt

class InterproError(Exception):
    pass

def _parse_none_if_minus(val):
    "parse fields that use '-' if there is no value"
    return None if val in ('-', '') else val.strip()

def _parse_status(val):
    "convert status to a boolean"
    return True if val == 'T' else False

# `TSV' does not contain a header
_tsvColumns = ("protein_accession", "sequence_digest", "sequence_length",
               "analysis", "signature_accession", "signature_description", "start", "stop",
               "score", "status", "date", "interpro_accession", "interpro_description",
               "go", "pathways")
_tsvTypeMap = {"sequence_length": int,
               "start": int,
               "stop": int,
               "score": _parse_none_if_minus,  # don't bother converting to float
               "status": _parse_status,
               "signature_description": _parse_none_if_minus,
               "interpro_accession": _parse_none_if_minus,
               "interpro_description": _parse_none_if_minus,
               "go": _parse_none_if_minus,
               "pathways": _parse_none_if_minus}

class InterproAnnot(TsvRow):
    """one interpro annotation record"""

    def short(self):
        desc = f"{self.protein_accession}/{self.analysis}:"
        if self.interpro_accession is not None:
            return f"{desc}: {self.interpro_accession}/{self.interpro_description}"
        else:
            return f"{desc}: {self.signature_description}"


class InterproAnnotTbl(list):
    """InterProScan analysis results from TSV output.

    Rows have the above columns from TSV, plus
    annotId - computed annotId with protein id and offset within proteins annotations
    """
    def __init__(self):
        self.byAcc = defaultdict(list)
        self.byAnnotId = {}

    def add(self, row):
        self.append(row)
        self.byAcc[row.protein_accession].append(row)
        row.annotId = annotIdFmt(row.protein_accession, len(self.byAcc[row.protein_accession]) - 1)
        assert row.annotId not in self.byAnnotId
        self.byAnnotId[row.annotId] = row

    def finish(self):
        self.byAcc.default_factory = None

    def getByAcc(self, acc):
        "Error if not found"
        try:
            return self.byAcc[acc]
        except Exception as ex:
            raise InterproError(f"InterPro acc not found {acc}") from ex

    def getByAnnotId(self, annotId):
        annot = self.byAnnotId.get(annotId)
        if annot is None:
            raise InterproError(f"InterPro annotId '{annotId}' not found in annotation table")
        return annot

def interproAnnotsLoad(interproTsv):
    """Load interproscan results TSV """

    interproTbl = InterproAnnotTbl()
    for row in TsvReader(interproTsv, columns=_tsvColumns, typeMap=_tsvTypeMap,
                         rowClass=InterproAnnot):
        interproTbl.add(row)
    interproTbl.finish()
    return interproTbl
