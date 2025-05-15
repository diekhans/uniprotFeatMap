from pycbio.tsv import TsvReader, strOrNoneType
from uniprotmap.geneset import GeneSet, geneSetLoadAnnotPsl, geneSetLoadAnnotGp

codingTranscriptTypes = frozenset(["protein_coding",
                                   "nonsense_mediated_decay",
                                   "non_stop_decay",
                                   "protein_coding_LoF",
                                   "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene",
                                   "IG_C_gene", "IG_J_gene", "IG_D_gene", "IG_V_gene"])

def _loadMetadata(geneSet, gencodeMetaTsv):
    """read from GENCODE metadata TSV, from UCSC GENCODE import"""
    for row in TsvReader(gencodeMetaTsv, typeMap={"proteinId": strOrNoneType}):
        if row.transcriptType in codingTranscriptTypes:
            geneSet.meta.addTranscript(row.geneId, row.geneName, row.geneType,
                                       row.transcriptId, row.transcriptType,
                                       row.proteinId)

def gencodeGeneSetFactory(geneSetName, *, geneSetMetadata=None, trans2GenomePsl=None, transGenomeGp=None, transFa=None):
    """load GENCODE data as requested"""

    geneSet = GeneSet(geneSetName)
    if geneSetMetadata is not None:
        _loadMetadata(geneSet, geneSetMetadata)
    if trans2GenomePsl is not None:
        geneSetLoadAnnotPsl(geneSet.data, trans2GenomePsl)
    if transGenomeGp is not None:
        geneSetLoadAnnotGp(geneSet.data, transGenomeGp)
    geneSet.transFa = transFa
    geneSet.finish()
    return geneSet
