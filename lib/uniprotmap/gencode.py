from pycbio.tsv import TsvReader, strOrNoneType
from uniprotmap.geneset import GeneSetMetadata, geneSetDataLoad

_codingTranscriptTypes = frozenset(["protein_coding",
                                    "nonsense_mediated_decay",
                                    "non_stop_decay",
                                    "protein_coding_LoF",
                                    "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene",
                                    "IG_C_gene", "IG_J_gene", "IG_D_gene", "IG_V_gene"])

def loadGencodeMetadata(gencodeMetaTsv):
    """read from GENCODE metadata TSV, from UCSC GENCODE import"""
    geneSetMeta = GeneSetMetadata()
    for row in TsvReader(gencodeMetaTsv, typeMap={"proteinId": strOrNoneType}):
        if row.transcriptType in _codingTranscriptTypes:
            geneSetMeta.addTranscript(row.geneId, row.geneName, row.geneType,
                                      row.transcriptId, row.transcriptType,
                                      row.proteinId)
    geneSetMeta.finish()
    return geneSetMeta

def loadGencodeData(gencodePsl, gencodeGp=None):
    """8read in GENCODE PSLs and option genePreds.  Should already be filtered
    for protein-coding"""
    return geneSetDataLoad(gencodePsl, gencodeGp)
