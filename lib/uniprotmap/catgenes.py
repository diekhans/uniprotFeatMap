"""
GeneSet implementation for CAT
"""
from uniprotmap.geneset import GeneSet, geneSetLoadAnnotPsl


def cat1GeneSetFactory(geneSetName, *, geneSetMetadata=None, transGenomePsl=None, transGenomeGp=None, transFa=None):
    """load CAT format 1 data as requested"""

    geneSet = GeneSet(geneSetName)
    if geneSetMetadata is not None:
        assert False, "cat1GeneSetFactory geneSetMetadata not implemented"
    if transGenomePsl is not None:
        geneSetLoadAnnotPsl(geneSet.data, transGenomePsl)
    if transGenomeGp is not None:
        assert False, "cat1GeneSetFactory transGenomeGp not implemented"
    geneSet.transFa = transFa
    geneSet.finish()
    return geneSet
