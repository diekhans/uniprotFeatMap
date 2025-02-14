"""
Common CLI parsing functions.
"""
from uniprotmap.geneset import GeneSetName

def cliAddGeneSetParameters(parser, *, inclMetadata=False, inclTransGenomePsl=False,
                            inclTransGenomeGp=False, inclTransFa=False):
    """add the requested position parameters to a parser. A geneSet name is
    always required."""
    parser.add_argument("geneSetName", type=GeneSetName, choices=GeneSetName,
                        help="""name of gene set""")
    if inclMetadata:
        parser.add_argument("geneSetMetadata",
                            help="""Metadata file for the GeneSet, format depends on which set (input)""")
    if inclTransGenomePsl:
        parser.add_argument("transGenomePsl",
                            help="""Transcript genome alignment; often from genePredToPsl. For cross-species mapping, this should be the target species transcript alignments. (input)""")
    if inclTransGenomeGp:
        parser.add_argument("transGenomeGp",
                            help="""Transcript genome alignment; often from genePredToPsl. For cross-species mapping, this should be the target species transcript alignments. (input)""")
    if inclTransFa:
        parser.add_argument("transFa",
                            help="""Transcript FASTA file. (input)""")
