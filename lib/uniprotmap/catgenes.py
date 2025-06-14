"""
GeneSet implementation for CAT
"""
from pycbio.tsv import TsvReader
from uniprotmap.geneset import GeneSet, geneSetLoadAnnotPsl
from uniprotmap.gencode import codingTranscriptTypes

def _loadMetadata(geneSet, geneSetMetadata):
    """read CAT metadata, which is part of the bigBed """
    for row in TsvReader(geneSetMetadata):
        if row.type in codingTranscriptTypes:
            geneSet.meta.addTranscript(row.geneName, row.name2, row.geneType,
                                       row.name, row.type,
                                       row.name + '_prot')

def cat1GeneSetFactory(geneSetName, *, geneSetMetadata=None, trans2GenomePslFile=None, transGenomeGpFile=None, transFa=None):
    """load CAT format 1 data as requested"""

    geneSet = GeneSet(geneSetName)
    if geneSetMetadata is not None:
        _loadMetadata(geneSet, geneSetMetadata)
    if trans2GenomePslFile is not None:
        geneSetLoadAnnotPsl(geneSet.data, trans2GenomePslFile)
    if transGenomeGpFile is not None:
        assert False, "cat1GeneSetFactory transGenomeGp not implemented"
    geneSet.transFa = transFa
    geneSet.finish()
    return geneSet
