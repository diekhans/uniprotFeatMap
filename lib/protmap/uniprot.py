"""
Reads files create by uniprotToTab
"""

import pandas as pd
from protmap import dropVersion

def splitMetaList(val):
    """ split strings like: ENST00000369413.8|ENST00000528909.1"""
    if val == "":
        return []
    return val.split('|')

def splitDropVersion(idents):
    return [dropVersion(id) for id in splitMetaList(idents)]

def dropUniportIsoformModifier(acc):
    return acc.split('-')[0]

class UniprotMetaTbl:
    """reads swissprot.9606.tab, trembl.9606.tab"""
    def __init__(self, uniprotMetaTsv):
        self.df = pd.read_table(uniprotMetaTsv, keep_default_na=False)

        # index by gene and transcript ids
        self.df['ensemblGeneAccs'] = self.df.ensemblGene.apply(splitDropVersion)
        self.df['ensemblTransAccs'] = self.df.ensemblTrans.apply(splitDropVersion)
        self.df['isoIds'] = self.df.isoIds.apply(splitMetaList)
        self.byGeneAccDf = self.df.explode('ensemblGeneAccs')
        self.byGeneAccDf.rename(columns={'ensemblGeneAccs': 'ensemblGeneAcc'}, inplace=True)
        self.byTranscriptAccDf = self.df.explode('ensemblTransAccs')
        self.byTranscriptAccDf.rename(columns={'ensemblTransAccs': 'ensemblTransAcc'}, inplace=True)

        # with only one isoform: acc == mainIsoAcc, and isoIds are empty
        # for multiple isoforms: mainIsoAcc is canonical and isoIds are other isoforms.
        # Build an index by isoId which is
        self.df['allIsoIds'] = self.df.apply(lambda row: row.isoIds + [row.mainIsoAcc], axis=1)
        self.byIsoId = self.df.explode('allIsoIds')
        self.byIsoId.rename(columns={'allIsoIds': 'isoId'}, inplace=True)

    def getByAcc(self, acc):
        return self.df[self.df.acc == acc].iloc[0]

    def getByIsoId(self, isoId):
        try:
            return self.byIsoId[self.byIsoId.isoId == isoId].iloc[0]
        except IndexError as ex:
            raise IndexError(f"isoId not found {isoId}") from ex

    def getTransMeta(self, transAcc):
        protMetas = self.byTranscriptAccDf[self.byTranscriptAccDf.ensemblTransAcc == transAcc]
        if len(protMetas) == 0:
            return None
        if len(protMetas) > 1:
            raise Exception(f"more than one protein for transcript {transAcc}, found: {protMetas.mainIsoAcc}")
        return protMetas.iloc[0]

    def getGeneAccMetas(self, geneAcc):
        protMetas = self.byGeneAccDf[self.byGeneAccDf.ensemblGeneAcc == geneAcc]
        if len(protMetas) == 0:
            return None
        return protMetas

    def getGeneNameMetas(self, geneName):
        protMetas = self.byGeneAccDf[self.byGeneAccDf.geneName == geneName]
        if len(protMetas) == 0:
            return None
        return protMetas


class UniprotAnnotTbl:
    """reads swissprot.9606.annots.tab, trembl.9606.annots.tab"""
    def __init__(self, uniprotAnnotsTsv):
        self.df = pd.read_table(uniprotAnnotsTsv, keep_default_na=False)


def canonicalAnnotEncode(canonId, annotIdx):
    return canonId + "#" + str(annotIdx)

def canonicalAnnotDecode(name):
    try:
        canonId, annotIdx = name.split('#')
        return canonId, int(annotIdx)
    except ValueError:
        raise ValueError(f"invalid encode canonId and annotIdx: '{name}'")
