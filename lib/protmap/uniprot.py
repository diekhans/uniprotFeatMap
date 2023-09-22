"""
Reads files create by uniprotToTab
"""

import pandas as pd
from protmap import dropVersion, buildDfUniqueIndex, buildDfMultiIndex

# WARNING: UniProt is 1-based, open-end


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
        if "orgName" not in self.df.columns:
            raise Exception(f"column 'orgName' not found, is this a UnIprot metadata TSV: {uniprotMetaTsv}")
        self.df.set_index('acc', inplace=True, drop=False, verify_integrity=True)
        self.geneNameIdx = buildDfMultiIndex(self.df, "geneName")

        # index by gene and transcript ids
        self.df['ensemblGeneAccs'] = self.df.ensemblGene.apply(splitDropVersion)
        self.df['ensemblTransAccs'] = self.df.ensemblTrans.apply(splitDropVersion)
        self.df['isoIds'] = self.df.isoIds.apply(splitMetaList)

        self.byGeneAccDf = self.df.explode('ensemblGeneAccs')
        self.byGeneAccDf.rename(columns={'ensemblGeneAccs': 'ensemblGeneAcc'}, inplace=True)
        self.byGeneAccDf.set_index('ensemblGeneAcc', inplace=True, drop=False, verify_integrity=True)

        self.byTranscriptAccDf = self.df.explode('ensemblTransAccs')
        self.byTranscriptAccDf.rename(columns={'ensemblTransAccs': 'ensemblTransAcc'}, inplace=True)
        self.byTranscriptAccDf.set_index('ensemblTransAcc', inplace=True, drop=False, verify_integrity=True)

        # with only one isoform: acc == mainIsoAcc, and isoIds are empty
        # for multiple isoforms: mainIsoAcc is canonical and isoIds are other isoforms.
        # Build an index by isoId which is
        self.df['allIsoIds'] = self.df.apply(lambda row: row.isoIds + [row.mainIsoAcc], axis=1)
        self.byIsoId = self.df.explode('allIsoIds')
        self.byIsoId.rename(columns={'allIsoIds': 'isoId'}, inplace=True)
        self.byIsoId.set_index('isoId', inplace=True, drop=False, verify_integrity=True)

    def getByAcc(self, acc):
        return self.df[self.df.index == acc].iloc[0]

    def getByIsoId(self, isoId):
        try:
            return self.byIsoId[self.byIsoId.index == isoId].iloc[0]
        except IndexError as ex:
            raise IndexError(f"isoId not found {isoId}") from ex

    def getTransMeta(self, transAcc):
        protMetas = self.byTranscriptAccDf[self.byTranscriptAccDf.index == transAcc]
        if len(protMetas) == 0:
            return None
        if len(protMetas) > 1:
            raise Exception(f"more than one protein for transcript {transAcc}, found: {protMetas.mainIsoAcc}")
        return protMetas.iloc[0]

    def getGeneAccMetas(self, geneAcc):
        protMetas = self.byGeneAccDf[self.byGeneAccDf.index == geneAcc]
        if len(protMetas) == 0:
            return None
        return protMetas

    def getGeneNameMetas(self, geneName):
        return self.geneNameIdx.get(geneName)


class UniprotAnnotTbl:
    """reads swissprot.9606.annots.tab, trembl.9606.annots.tab"""
    def __init__(self, uniprotAnnotsTsv):
        self.df = pd.read_table(uniprotAnnotsTsv, keep_default_na=False)
        if "featType" not in self.df.columns:
            raise Exception("column 'featType' not found, is this a UnIprot annotations TSV: " + uniprotAnnotsTsv)

        # add a unique, reproducible annotation id in the form mainIsoAcc#annotIdx, where the annotIdx is
        # the relative index of the annotation for that acc.
        self.df["annotId"] = self.df.mainIsoAcc.astype(str) + "#" + self.df.groupby("mainIsoAcc").cumcount().astype(str)
        #self.df.set_index('annotId', inplace=True, drop=False, verify_integrity=True)
        self.annotIdIdx = buildDfUniqueIndex(self.df, "annotId")

    def getByAnnotId(self, annotId):
        #annots = self.df[self.df.index == annotId]
        #if len(annots) == 0:
        #    raise Exception(f"annotId '{annotId}' not found in annotation table")
        #return annots.iloc[0]

        annot = self.annotIdIdx.get(annotId)
        if annot is None:
            raise Exception(f"annotId '{annotId}' not found in annotation table")
        return annot
